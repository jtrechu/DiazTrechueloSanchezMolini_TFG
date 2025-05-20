# Load the data used for modeling from an external R script
source("TableData.R")
samples <- tableData()  # Generate the dataset

# Set Elastic Net mixing parameter
alpha = 0.5

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
set.seed(123)

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

# Create 5-fold cross-validation groups based on the 'Group' variable
folds <- createFolds(samples$Group, k = 5, list = TRUE, returnTrain = FALSE)

# Define training set for fold g
train_indices <- setdiff(1:nrow(samples), folds[[g]])

# Initialize matrix to store MSEs for different gamma values and folds
MSE.df <- matrix(numeric(0), nrow = 4 , ncol = 5)
betas_fold <- list()

# Cross-validation loop over 5 folds
for (g in 1:5){
  betas <- matrix(numeric(0), nrow = 21, ncol = 4)  # Store betas per gamma
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
  lambda <- vector_lambdas[g]
  
  # Fit an Elastic Net model with alpha and lambda
  en <- glmnet(Xtrain[,-1], ytrain, alpha = alpha, lambda = lambda / 2, standardize = FALSE)
  en <- c(en$a0, as.matrix(en$beta))  # Combine intercept and coefficients
  
  # Define gamma values for cost-sensitive constraint
  gammas <- c(0.03, 0.05, 0.1, 0.15)
  MSE <- c()
  
  for (gamma in gammas){
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
    
    # Evaluate test error (MSE) for this model
    Xtest <- samples[folds[[g]],-c(22,23)]
    ytest <- samples[folds[[g]],22]
    residuos <- as.matrix(Xtest) %*% betas_aux - ytest
    MSE_aux <- sum(residuos^2) / length(residuos)
    MSE <- c(MSE, MSE_aux)
    
    # Store coefficients for this gamma value
    betas[,which(gammas == gamma)] <- betas_aux
  }
  
  # Store MSE values and coefficients for this fold
  MSE.df[,g] <- MSE
  betas_fold[[g]] <- betas
}

# Identify best performing gamma (one with minimum MSE) for each gamma index
cvIndices <- apply(MSE.df, 1, which.min)
cvIndices

# ----- Re-Run previous code for Final model fitting 
# ----- to obtain MSE for different groups and parameters
length(which(beta>1e-7))
residuos <- as.matrix(Xtrain)%*%beta - ytrain
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_1)%*%beta - y_1
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_2)%*%beta - y_2
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_3)%*%beta - y_3
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_4)%*%beta - y_4
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_5)%*%beta - y_5
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_6)%*%beta - y_6
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_7)%*%beta - y_7
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_8)%*%beta - y_8
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
