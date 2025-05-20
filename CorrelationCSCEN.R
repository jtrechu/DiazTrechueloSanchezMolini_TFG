options(digits = 11)

# Load required libraries
library(gurobi)
library(ggplot2)
library(latex2exp)
library(tidyr)
library(dplyr)

# Load user-defined functions
source('ECM_star_function.R')
source('tau_min_function_2.R')
source('tau_max_function.R')
source('CorrelationDataSet.R')

run_constrained_lasso <- function(alpha, p = 0.99, ext_lambdas = 30, grid_lambda = 0.1, grid_tau = 0.1, taus = c(50)) {
  vector_lambdas <- seq(0, ext_lambdas, by = grid_lambda)
  
  # Generate dataset
  df <- generate_dataset(p, 100, 80)
  df <- cbind(Intercept = 1, df)
  
  # Extract design matrix and response
  X <- as.matrix(df[, -c(2, 13)])
  y <- df[, 2]
  
  # Subset for group 2
  X_2 <- as.matrix(df[df[, 13] == 2, -c(2, 13)])
  y_2 <- df[df[, 13] == 2, 2]
  
  # Compute baseline MSE for constraints
  constraints_datasets <- list(as.matrix(cbind(y_2, X_2)))
  f0_aux <- ECM_star_function(X_2, y_2)
  ECM_star <- as.numeric(f0_aux[2])
  
  # Optimization output containers
  betas <- list()
  valorobjetivo <- list()
  status_vector <- list()
  j <- 0
  
  for (lambda in vector_lambdas) {
    j <- j + 1
    betas_aux <- c()
    valorobjetivo_aux <- c()
    status_vector_aux <- c()
    
    for (tau in taus) {
      # Create penalty matrix A (for elastic net penalty)
      A_aux <- cbind(rep(0, ncol(X) - 1), diag(ncol(X) - 1))
      A <- cbind(A_aux, -diag(nrow(A_aux)), diag(nrow(A_aux)))
      
      # Objective function components
      objcon <- (1 / nrow(X)) * t(y) %*% y
      Q <- (1 / nrow(X)) * rbind(
        cbind(t(X) %*% X + lambda * (1 - alpha) * (t(A_aux) %*% A_aux),
              matrix(0, ncol = 2 * nrow(A), nrow = ncol(X))),
        matrix(0, nrow = 2 * nrow(A), ncol = ncol(X) + 2 * nrow(A))
      )
      obj <- c(rep(0, ncol(X)), rep(alpha * lambda, 2 * nrow(A))) -
        (2 / nrow(X)) * t(y) %*% cbind(X, matrix(0, nrow = nrow(X), ncol = 2 * nrow(A)))
      
      # Quadratic constraints
      quadcon <- list()
      for (i in seq_along(constraints_datasets)) {
        data_i <- constraints_datasets[[i]]
        X_i <- data_i[, -1]
        y_i <- data_i[, 1]
        Qc_i <- (1 / nrow(X_i)) * rbind(
          cbind(t(X_i) %*% X_i, matrix(0, ncol = 2 * nrow(A), nrow = ncol(X_i))),
          matrix(0, ncol = ncol(X_i) + 2 * nrow(A), nrow = 2 * nrow(A))
        )
        q_i <- -(2 / nrow(X_i)) * t(y_i) %*% cbind(X_i, matrix(0, nrow = nrow(X_i), ncol = 2 * nrow(A)))
        rhs_i <- (1 + tau) * ECM_star - (1 / nrow(X_i)) * t(y_i) %*% y_i
        quadcon[[i]] <- list(Qc = Qc_i, q = as.numeric(q_i), rhs = rhs_i, sense = "<")
      }
      
      # Define model for Gurobi
      model <- list(
        A = A, rhs = 0, sense = "=", obj = obj, Q = Q,
        vtype = "C", modelsense = "min", lb = c(rep(-Inf, ncol(X)), rep(0, 2 * nrow(A))),
        objcon = objcon, quadcon = quadcon
      )
      
      # Run optimization
      result <- gurobi(model)
      status_vector_aux <- c(status_vector_aux, result$status)
      valorobjetivo_aux <- c(valorobjetivo_aux, result$objval)
      betas_aux <- rbind(betas_aux, result$x[1:ncol(X)])
    }
    
    betas[[j]] <- betas_aux
    valorobjetivo[[j]] <- valorobjetivo_aux
    status_vector[[j]] <- status_vector_aux
  }
  
  list(betas = betas, valorobjetivo = valorobjetivo, status = status_vector)
}

# Run the model for alpha = 0.25
result_025 <- run_constrained_lasso(alpha = 0.25)

# Process results
lambda <- seq(0, 30, by = 0.1)
beta_matrix_025 <- do.call(rbind, result_025$betas)
beta_matrix_025 <- cbind(lambda, beta_matrix_025[, c(2, 3, 4)])
colnames(beta_matrix_025) <- c("lambda", paste0("beta", 1:3))
df_025 <- as.data.frame(beta_matrix_025)
df_025$DifEN <- abs(df_025$beta1 - df_025$beta2)

# Run for alpha = 0.5
result_050 <- run_constrained_lasso(alpha = 0.5)
beta_matrix_050 <- do.call(rbind, result_050$betas)
beta_matrix_050 <- cbind(lambda, beta_matrix_050[, c(2, 3, 4)])
colnames(beta_matrix_050) <- c("lambda", paste0("beta", 1:3))
df_050 <- as.data.frame(beta_matrix_050)
df_050$DifEN <- abs(df_050$beta1 - df_050$beta2)

# Run for alpha = 0.75
result_075 <- run_constrained_lasso(alpha = 0.75)
beta_matrix_075 <- do.call(rbind, result_075$betas)
beta_matrix_075 <- cbind(lambda, beta_matrix_075[, c(2, 3, 4)])
colnames(beta_matrix_075) <- c("lambda", paste0("beta", 1:3))
df_075 <- as.data.frame(beta_matrix_075)
df_075$DifEN <- abs(df_075$beta1 - df_075$beta2)

# Run for alpha = 0.75
result_1 <- run_constrained_lasso(alpha = 1)
beta_matrix_1 <- do.call(rbind, result_1$betas)
beta_matrix_1 <- cbind(lambda, beta_matrix_1[, c(2, 3, 4)])
colnames(beta_matrix_1) <- c("lambda", paste0("beta", 1:3))
df_1 <- as.data.frame(beta_matrix_1)
df_1$DifEN <- abs(df_1$beta1 - df_1$beta2)

library(ggplot2)
library(dplyr)
library(tidyr)
library(latex2exp)

# Combine all data frames with an alpha identifier
df_025$alpha <- 0.25
df_050$alpha <- 0.5
df_075$alpha <- 0.75
df_1$alpha <- 1

# Combine into one data frame
df_all <- bind_rows(df_025, df_050, df_075,df_1)

# Plot |β1 - β2| over lambda for all alpha values
ggplot(df_all, aes(x = lambda, y = DifEN, color = factor(alpha))) +
  geom_line(size = 1) +
  labs(
    title = TeX("Difference $|\\beta_1 - \\beta_2|$ vs. $\\lambda$ for different $\\alpha$"),
    x = TeX("$\\lambda$"),
    y = TeX("$|\\beta_1 - \\beta_2|$"),
    color = TeX("$\\alpha$")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )



