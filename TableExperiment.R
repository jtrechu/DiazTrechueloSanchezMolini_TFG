# Load the data used for modeling from an external R script
source("TableData.R") # Generate the dataset
resultados_totales <- data.frame()
# Set Elastic Net mixing parameter
for (alpha in c(1,0.75,0.5,0.25)){

# Configure number display options
options(digits = 11)

# Load required libraries
library(gurobi)
path_directory = "C:/Users/jtrec/OneDrive/Desktop/Codigo Jaime/"

# Define the alpha value(s) for EN
vector_alphas <- c(alpha)

# Maximum lambda value for EN regularization
ext_lambdas <- 12

# Initialize list to store coefficients from each fold
betas_fold <- list()

# Create a sequence of lambda values for testing, dropping the first one (zero)
vector_lambdas <- seq(0, ext_lambdas, length.out = 6)[-1]

# Load additional libraries for modeling and cross-validation
library(genridge)
library(glmnet)
library(caret)

# Set random seed for reproducibility
set.seed(2025)
samples <- tableData() 
# Prepare feature matrix and response variable
X <- samples[,-c(22,23)]  # Drop response and group columns
X <- as.matrix(X)
y <- samples[,22]         # Response variable

# Split data into subgroups (1 to 8) based on the 'Group' variable
X_1 <- samples[which(samples[,"Group"]==1),-c(22,23)]
y_1 <- samples[which(samples[,"Group"]==1),22]
X_2 <- samples[which(samples[,"Group"]==2),-c(22,23)]
y_2 <- samples[which(samples[,"Group"]==2),22]
X_3 <- samples[which(samples[,"Group"]==3),-c(22,23)]
y_3 <- samples[which(samples[,"Group"]==3),22]
X_4 <- samples[which(samples[,"Group"]==4),-c(22,23)]
y_4 <- samples[which(samples[,"Group"]==4),22]
X_5 <- samples[which(samples[,"Group"]==5),-c(22,23)]
y_5 <- samples[which(samples[,"Group"]==5),22]
X_6 <- samples[which(samples[,"Group"]==6),-c(22,23)]
y_6 <- samples[which(samples[,"Group"]==6),22]
X_7 <- samples[which(samples[,"Group"]==7),-c(22,23)]
y_7 <- samples[which(samples[,"Group"]==7),22]
X_8 <- samples[which(samples[,"Group"]==8),-c(22,23)]
y_8 <- samples[which(samples[,"Group"]==8),22]

# Create 2-folds train and test
library(caret)
train_indices <- createDataPartition(samples$Group, p = 0.75, list = FALSE)

# Initialize matrix to store MSEs for different gamma values and folds
MSE.df <- matrix(numeric(0), nrow = 4 , ncol = 5)
betas_fold <- list()

# Cross-validation loop over 5 folds

  betas <- matrix(numeric(0), nrow = 21, ncol = 5)  # Store betas per gamma
  Xtrain <- samples[train_indices,-c(22,23)]
  Xtrain <- as.matrix(Xtrain)
  ytrain <- samples[train_indices,22]
  
  # Prepare group-specific data from the training set
  X_1train <- samples[which(samples[train_indices,"Group"]==1),-c(22,23)]
  y_1train <- samples[which(samples[train_indices,"Group"]==1),22]
  X_2train <- samples[which(samples[train_indices,"Group"]==2),-c(22,23)]
  y_2train <- samples[which(samples[train_indices,"Group"]==2),22]
  X_3train <- samples[which(samples[train_indices,"Group"]==3),-c(22,23)]
  y_3train <- samples[which(samples[train_indices,"Group"]==3),22]
  
  # Select the lambda for this fold
  lambda <- 2.4
  
  # Fit an Elastic Net model with alpha and lambda
  en <- glmnet(Xtrain[,-1], ytrain, alpha = alpha, lambda = lambda / 2, standardize = FALSE)
  en <- c(en$a0, as.matrix(en$beta))  # Combine intercept and coefficients
  
  # Define gamma values for cost-sensitive constraint
  gammas <- c(0.03, 0.05, 0.1, 0.15)
  MSE <- c()
  
  for (k in (1:length(gammas))){
    gamma = gammas[k]
    # Compute adjusted MSE terms for group 1 to 3
    f1 <- (1-gamma)*mean((y_1train - as.matrix(X_1train)%*%en)^2)
    f2 <- (1-gamma)*mean((y_2train - as.matrix(X_2train)%*%en)^2)
    f3 <- (1-gamma)*mean((y_3train - as.matrix(X_3train)%*%en)^2)
    f <- c(f1,f2,f3)
    
    # Optimization setup
    j = 0
    valorobjetivo <- list()
    status_vector <- list()
    j = j + 1
    valorobjetivo_aux <- c()
    betas_aux <- c()
    status_vector_aux <- c()
    
    # Identity structure for EN penalty
    A_aux <- cbind(rep(0,ncol(Xtrain)-1), diag(ncol(Xtrain)-1))
    A <- as.matrix(cbind(A_aux, -diag(nrow(A_aux)), diag(nrow(A_aux))))
    
    # Define constraint direction, objective, and quadratic penalty
    rhs <- 0
    sense <- '='
    objcon <- (1/nrow(Xtrain)) * t(ytrain) %*% ytrain
    
    # Construct the full quadratic matrix for optimization
    matrizcompuesta <- (1/nrow(Xtrain)) * as.matrix(rbind(
      cbind(t(Xtrain)%*%Xtrain + lambda*(1-alpha)*(t(A_aux)%*%A_aux), matrix(rep(0,ncol(Xtrain)*2*nrow(A)), ncol=2*nrow(A))),
      matrix(rep(0, 2*nrow(A)*(ncol(Xtrain)+2*nrow(A))), nrow=2*nrow(A))
    ))
    
    # Linear part of the objective function
    obj <- c(rep(0, ncol(Xtrain)), rep(alpha*lambda, 2*nrow(A))) - (2/nrow(Xtrain)) * t(ytrain) %*% cbind(Xtrain, matrix(rep(0, 2*nrow(A)*nrow(Xtrain)), ncol=2*nrow(A)))
    
    Q <- matrizcompuesta
    vtype <- 'C'
    modelsense <- 'min'
    modelname <- 'CLassoMatrizId'
    lb <- c(rep(-Inf,ncol(Xtrain)), rep(0,2*nrow(A)))
    
    # Define empty list to store quadratic constraints
    quadcon <- list()
    Cuadratica <- list()
    
    # Build constraint datasets for groups 1, 2, and 3
    constraints_datasets <- list(list())
    constraints_datasets[[1]] <- as.matrix(cbind(y_1train,X_1train))
    constraints_datasets[[2]] <- as.matrix(cbind(y_2train,X_2train))
    constraints_datasets[[3]] <- as.matrix(cbind(y_3train,X_3train))
    constraint_number <- length(constraints_datasets)
    
    # Define the quadratic constraints (cost-sensitive constraints)
    for (i in 1:constraint_number){
      matrizcompuestarestriccion <- as.matrix(rbind(
        cbind(t(constraints_datasets[[i]][,-1])%*%constraints_datasets[[i]][,-1], matrix(rep(0, ncol(constraints_datasets[[i]][,-1])*2*nrow(A)), ncol=2*nrow(A))),
        matrix(rep(0, 2*nrow(A)*(ncol(constraints_datasets[[i]][,-1])+2*nrow(A))), nrow=2*nrow(A))
      ))
      Cuadratica$Qc[[i]] = (1/nrow(constraints_datasets[[i]][,-1])) * matrizcompuestarestriccion
      Cuadratica$q[[i]] = -(2/nrow(constraints_datasets[[i]][,-1])) * t(constraints_datasets[[i]][,1]) %*% cbind(constraints_datasets[[i]][,-1], matrix(rep(0, 2*nrow(A)*nrow(constraints_datasets[[i]][,-1])), ncol=2*nrow(A)))
      Cuadratica$rhs[[i]] = f[i] - (1/nrow(constraints_datasets[[i]][,-1])) * t(constraints_datasets[[i]][,1]) %*% constraints_datasets[[i]][,1]
      Cuadratica$sense[[i]] <- '<'
      
      quadcon[[i]] <- list(Qc = Cuadratica$Qc[[i]], q = as.numeric(Cuadratica$q[[i]]), rhs = Cuadratica$rhs[[i]], sense = Cuadratica$sense[[i]])
    }
    
    # Build Gurobi model
    model <- list(A=A, rhs=rhs, sense=sense, obj=obj, Q=Q, vtype=vtype, modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
    
    # Solve the constrained optimization
    result <- gurobi(model)
    
    # Store solution status and objective value
    status_vector_aux <- c(status_vector_aux, result$status)
    valorobjetivo_aux1 <- result$objval
    
    # Extract beta coefficients from solution
    betas_aux <- result$x[1:ncol(Xtrain)]
    # Store coefficients for this gamma value
    betas[,k+1] <- betas_aux
  }
  betas[,1] <- as.numeric(en)

# ----- Re-Run previous code for Final model fitting 
# ----- to obtain MSE for different groups and parameters

Xtest <- X[-train_indices,]
ytest <- y[-train_indices]

Xtest_1 <- X_1[-train_indices,]
ytest_1 <- y_1[-train_indices]

Xtest_2 <- X_2[-train_indices, ]
ytest_2 <- y_2[-train_indices]

Xtest_3 <- X_3[-train_indices, ]
ytest_3 <- y_3[-train_indices]

Xtest_4 <- X_4[-train_indices, ]
ytest_4 <- y_4[-train_indices]

Xtest_5 <- X_5[-train_indices, ]
ytest_5 <- y_5[-train_indices]

Xtest_6 <- X_6[-train_indices, ]
ytest_6 <- y_6[-train_indices]

Xtest_7 <- X_7[-train_indices, ]
ytest_7 <- y_7[-train_indices]

Xtest_8 <- X_8[-train_indices, ]
ytest_8 <- y_8[-train_indices]

# Vector de valores de gamma
gammas <- c(0, 0.03, 0.05, 0.1, 0.15)

# Inicializar listas para almacenar resultados
mse_matrix <- matrix(NA, nrow = length(gammas), ncol = 8)
mse_global <- numeric(length(gammas))  # <- NUEVO
num_coef_vec <- numeric(length(gammas))

# Rellenar matrices y vectores
for (i in 1:length(gammas)) {
  beta <- betas[, i]
  num_coef <- length(which(beta > 1e-6))
  num_coef_vec[i] <- num_coef
  
  # Calcular MSE global
  residuos_global <- as.matrix(Xtest) %*% beta - ytest
  mse_global[i] <- sum(residuos_global^2) / length(residuos_global)
  
  for (j in 1:8) {
    Xtest_j <- get(paste0("Xtest_", j))
    ytest_j <- get(paste0("ytest_", j))
    residuos <- as.matrix(Xtest_j) %*% beta - ytest_j
    mse_matrix[i, j] <- sum(residuos^2) / length(residuos)
  }
}

# Crear el data.frame final con MSE global al principio
colnames(mse_matrix) <- paste0("MSE_test_", 1:8)
resultados <- data.frame(
  alpha = alpha,
  Gamma = gammas,
  NumCoeficientes = num_coef_vec,
  `MSE global` = mse_global,   # <- Aquí se añade
  mse_matrix
)

# Acumular resultados por alpha
resultados_totales <- rbind(resultados_totales, resultados)}
# Renombrar columnas para presentación LaTeX si se desea
colnames(resultados_totales) <- c(
  "$\\alpha$", "$\\gamma$","Vars","MSE",
  paste0("MSE ", 1:8)
)
# Crear tabla LaTeX
kable(resultados_totales, format = "latex", booktabs = TRUE, digits = 5,
      caption = "MSE por grupo para distintos valores de $\\alpha$ y $\\gamma$") %>%
  kable_styling(latex_options = c()) %>%
  collapse_rows(columns = 1, valign = "middle")  # <- aquí agrupamos alpha cada 4 filas
