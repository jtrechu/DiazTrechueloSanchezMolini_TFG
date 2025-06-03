# Load dataset-generating script
source("TableData.R")
resultados_totales <- data.frame()

# Libraries
library(gurobi)
library(genridge)
library(glmnet)
library(caret)

# Elastic Net mixing parameters
for (alpha in c(1, 0.75, 0.5, 0.25)) {
  set.seed(2025)
  options(digits = 11)
  
  samples <- tableData()
  
  # Define feature matrix X and response y
  X <- as.matrix(samples[,-c(22,23)])
  y <- samples[,22]
  
  # Training indices (75%)
  train_indices <- createDataPartition(samples$Group, p = 0.75, list = FALSE)
  
  Xtrain <- as.matrix(samples[train_indices, -c(22,23)])
  ytrain <- samples[train_indices, 22]
  lambda <- 2.4
  
  # Fit unconstrained Elastic Net
  en_model <- glmnet(Xtrain[,-1], ytrain, alpha = alpha, lambda = lambda / 2, standardize = FALSE)
  en_coef <- c(en_model$a0, as.matrix(en_model$beta))
  
  # Gamma values for fairness constraints
  gammas <- c(0, 0.03, 0.05, 0.1, 0.15)
  betas <- matrix(0, nrow = 21, ncol = length(gammas))
  betas[,1] <- as.numeric(en_coef)  # First column is baseline EN
  
  # Prepare group-specific training data
  for (k in 2:length(gammas)) {
    gamma <- gammas[k]
    
    group_datasets <- lapply(1:3, function(g) {
      idx <- which(samples$Group[train_indices] == g)
      cbind(samples[train_indices[idx], 22], as.matrix(samples[train_indices[idx], -c(22,23)]))
    })
    
    # Compute unconstrained group losses
    f <- sapply(group_datasets, function(d) {
      (1 - gamma) * mean((d[,1] - as.matrix(d[,-1]) %*% en_coef)^2)
    })
    
    # Build EN penalty structure
    A_aux <- cbind(rep(0, ncol(Xtrain)-1), diag(ncol(Xtrain)-1))
    A <- cbind(A_aux, -diag(nrow(A_aux)), diag(nrow(A_aux)))
    rhs <- 0
    sense <- '='
    objcon <- (1/nrow(Xtrain)) * sum(ytrain^2)
    
    Q <- (1/nrow(Xtrain)) * rbind(
      cbind(t(Xtrain) %*% Xtrain + lambda * (1 - alpha) * t(A_aux) %*% A_aux,
            matrix(0, ncol = 2 * nrow(A), nrow = ncol(Xtrain))),
      matrix(0, ncol = ncol(Xtrain) + 2 * nrow(A), nrow = 2 * nrow(A))
    )
    
    obj <- c(rep(0, ncol(Xtrain)), rep(alpha * lambda, 2 * nrow(A))) -
      (2 / nrow(Xtrain)) * t(ytrain) %*% cbind(Xtrain, matrix(0, nrow(Xtrain), 2 * nrow(A)))
    
    lb <- c(rep(-Inf, ncol(Xtrain)), rep(0, 2 * nrow(A)))
    
    # Quadratic constraints per group
    quadcon <- lapply(1:3, function(i) {
      d <- group_datasets[[i]][,-1]
      n <- nrow(d)
      Qc <- (1/n) * rbind(
        cbind(t(d) %*% d, matrix(0, ncol = 2 * nrow(A), nrow = ncol(d))),
        matrix(0, nrow = 2 * nrow(A), ncol = ncol(d) + 2 * nrow(A))
      )
      q <- -(2/n) * t(group_datasets[[i]][,1]) %*% cbind(d, matrix(0, nrow(d), 2 * nrow(A)))
      rhs <- f[i] - (1/n) * sum(group_datasets[[i]][,1]^2)
      list(Qc = Qc, q = as.numeric(q), rhs = rhs, sense = '<')
    })
    
    model <- list(
      A = A, rhs = rhs, sense = sense,
      obj = obj, Q = Q, vtype = 'C', modelsense = 'min', lb = lb,
      objcon = objcon, quadcon = quadcon
    )
    
    result <- gurobi(model)
    betas[,k] <- result$x[1:ncol(Xtrain)]
  }
  
  # Prepare test sets
  Xtest <- X[-train_indices,]
  ytest <- y[-train_indices]
  mse_matrix <- matrix(NA, nrow = length(gammas), ncol = 8)
  mse_global <- numeric(length(gammas))
  num_coef_vec <- numeric(length(gammas))
  
  for (i in seq_along(gammas)) {
    beta <- betas[,i]
    preds <- Xtest %*% beta
    mse_global[i] <- mean((ytest - preds)^2)
    num_coef_vec[i] <- sum(abs(beta) > 1e-6)
    
    for (j in 1:8) {
      group_idx <- which(samples$Group[-train_indices] == j)
      X_j <- Xtest[group_idx, ]
      y_j <- ytest[group_idx]
      mse_matrix[i,j] <- mean((X_j %*% beta - y_j)^2)
    }
  }
  
  resultados <- data.frame(
    alpha = alpha,
    Gamma = gammas,
    Num_Coefficients = num_coef_vec,
    Global_MSE = mse_global,
    mse_matrix
  )
  
  resultados_totales <- rbind(resultados_totales, resultados)
}

colnames(resultados_totales) <- c(
  "$\\alpha$", "$\\gamma$", "Vars", "MSE",
  paste0("MSE ", 1:8)
)

# LaTeX table generation
library(knitr)
library(kableExtra)
kable(resultados_totales, format = "latex", booktabs = TRUE, digits = 5,
      caption = "MSE by group for different values of $\\alpha$ and $\\gamma$") %>%
  kable_styling(latex_options = c()) %>%
  collapse_rows(columns = 1, valign = "middle")
