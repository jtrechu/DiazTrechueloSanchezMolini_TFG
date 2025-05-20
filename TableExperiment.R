source("TableData.R")
samples <- tableData()
alpha=0.5
options(digits = 11)
library(gurobi)
path_directory = "C:/Users/jtrec/OneDrive/Desktop/Codigo Jaime/"
vector_alphas <- c(alpha)
ext_lambdas <- 12
betas_fold <- list()
vector_lambdas <- seq(0,ext_lambdas,length.out=6)[-1]
library(gurobi)
library(genridge)
library(glmnet)
library(caret)
set.seed(123)

X <- samples[,-c(22,23)]
X <- as.matrix(X)
y <-samples[,22]
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
folds <- createFolds(samples$Group, k = 5, list = TRUE, returnTrain = FALSE)
train_indices <- setdiff(1:nrow(samples), folds[[g]])
MSE.df <- matrix(numeric(0), nrow = 4 , ncol = 5)
betas_fold <- list()
for (g in 1:5){
  betas <- matrix(numeric(0), nrow = 21, ncol = 4)
  Xtrain <- samples[train_indices,-c(22,23)]
  Xtrain <- as.matrix(Xtrain)
  ytrain <-samples[train_indices,22]
  X_1train <- (samples[which(samples[train_indices,"Group"]==1),-c(22,23)])
  y_1train <- samples[which(samples[train_indices,"Group"]==1),22]
  X_2train <- samples[which(samples[train_indices,"Group"]==2),-c(22,23)]
  y_2train <- samples[which(samples[train_indices,"Group"]==2),22]
  X_3train <- samples[which(samples[train_indices,"Group"]==3),-c(22,23)]
  y_3train<- samples[which(samples[train_indices,"Group"]==3),22]
  lambda<-vector_lambdas[g]
  en<-glmnet(Xtrain[,-1],ytrain,alpha = alpha,lambda = vector_lambdas[g]/2,standardize = F)
  en <- c(en$a0, as.matrix(en$beta))
  gammas <- c(0.03,0.05,0.1,0.15)
  MSE <- c()
  for (gamma in gammas ){
    
    f1 <- (1-gamma)*mean((y_1train - as.matrix(X_1train)%*%en)^2)
    f2 <- (1-gamma)*mean((y_2train - as.matrix(X_2train)%*%en)^2)
    f3 <- (1-gamma)*mean((y_3train - as.matrix(X_3train)%*%en)^2)
    f <- c(f1,f2,f3)
    j=0
    valorobjetivo <- list()
    status_vector <- list()
    j=j+1
    valorobjetivo_aux <- c()
    betas_aux <-c()
    status_vector_aux <- c()
    A_aux <- cbind(rep(0,ncol(Xtrain)-1), diag(ncol(Xtrain)-1)) #Tendra que cambiarse segun el tipo de problema
    A <- as.matrix(cbind(A_aux, -diag(nrow(A_aux)), diag(nrow(A_aux))))
    rhs        <- 0
    sense      <- '='
    objcon <- (1/nrow(Xtrain))*t(ytrain)%*%ytrain
    matrizcompuesta <- (1/nrow(Xtrain))*as.matrix(rbind(cbind(t(Xtrain)%*%Xtrain+lambda*(1-alpha)*(t(A_aux)%*%A_aux),matrix(rep(0,ncol(Xtrain)*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(Xtrain)+2*nrow(A))),nrow=2*nrow(A))))
    obj <- c(rep(0, ncol(Xtrain)), rep(alpha*lambda, 2*nrow(A))) - (2/nrow(Xtrain))*t(ytrain)%*%cbind(Xtrain,matrix(rep(0, 2*nrow(A)*nrow(Xtrain)),ncol=2*nrow(A)))
    Q <- matrizcompuesta
    vtype <- 'C'
    modelsense <- 'min'
    modelname <- 'CLassoMatrizId'
    lb <- c(rep(-Inf,ncol(Xtrain)), rep(0,2*nrow(A)))
    quadcon <- list()
    Cuadratica <- list()
    constraints_datasets <- list(list())
    constraints_datasets[[1]] <- cbind(y_1train,X_1train)
    constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
    constraints_datasets[[2]] <- cbind(y_2train,X_2train)
    constraints_datasets[[2]] <- as.matrix(constraints_datasets[[2]])
    constraints_datasets[[3]] <- cbind(y_3train,X_3train)
    constraints_datasets[[3]] <- as.matrix(constraints_datasets[[3]])
    constraint_number <- length(constraints_datasets)
    
    for (i in 1:constraint_number){
      matrizcompuestarestriccion <- as.matrix(rbind(cbind(t(constraints_datasets[[i]][,-1])%*%constraints_datasets[[i]][,-1],matrix(rep(0,ncol(constraints_datasets[[i]][,-1])*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(constraints_datasets[[i]][,-1])+2*nrow(A))),nrow=2*nrow(A))))
      Cuadratica$Qc[[i]]= (1/nrow(constraints_datasets[[i]][,-1]))*matrizcompuestarestriccion
      Cuadratica$q[[i]] = -(2/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%cbind(constraints_datasets[[i]][,-1],matrix(rep(0, 2*nrow(A)*nrow(constraints_datasets[[i]][,-1])),ncol=2*nrow(A)))
      Cuadratica$rhs[[i]] =  f[i] - (1/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%constraints_datasets[[i]][,1]
      Cuadratica$sense[[i]] <- '<'
      quadcon[[i]] <- list(Qc =  Cuadratica$Qc[[i]], q =  as.numeric(Cuadratica$q[[i]]), rhs= Cuadratica$rhs[[i]], sense = Cuadratica$sense[[i]])}
    
    model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
    result <- gurobi(model)
    status_vector_aux <- c(status_vector_aux, result$status)
    valorobjetivo_aux1 <- result$objval
    betas_aux <- result$x[1:ncol(Xtrain)]
    Xtest <- samples[folds[[g]],-c(22,23)]
    ytest <-samples[folds[[g]],22]
    residuos <- as.matrix(Xtest)%*%betas_aux - ytest
    MSE_aux <- sum(residuos^2)/length(residuos)
    MSE <- c(MSE,MSE_aux)
    betas[,which(gammas == gamma)] <- betas_aux
  }
  MSE.df[,g] <- MSE
  betas_fold[[g]] <- betas
}


cvIndices <- apply(MSE.df, 1, which.min)   
cvIndices 

##OBTAIN THE CV PARAMETERS

betas <- matrix(numeric(0), nrow = 21, ncol = 4)
Xtrain <- samples[,-c(22,23)]
Xtrain <- as.matrix(Xtrain)
ytrain <-samples[,22]
X_1train <- (samples[which(samples[,"Group"]==1),-c(22,23)])
y_1train <- samples[which(samples[,"Group"]==1),22]
X_2train <- samples[which(samples[,"Group"]==2),-c(22,23)]
y_2train <- samples[which(samples[,"Group"]==2),22]
X_3train <- samples[which(samples[,"Group"]==3),-c(22,23)]
y_3train<- samples[which(samples[,"Group"]==3),22]
lambda<-vector_lambdas[1]
en<-glmnet(Xtrain[,-1],ytrain,alpha = alpha,lambda = vector_lambdas[1]/2,standardize = F)
en <- c(en$a0, as.matrix(en$beta))
gammas <- c(0.03,0.05,0.1,0.15)
MSE <- c()
for (gamma in gammas ){
  
  f1 <- (1-gamma)*mean((y_1train - as.matrix(X_1train)%*%en)^2)
  f2 <- (1-gamma)*mean((y_2train - as.matrix(X_2train)%*%en)^2)
  f3 <- (1-gamma)*mean((y_3train - as.matrix(X_3train)%*%en)^2)
  f <- c(f1,f2,f3)
  j=0
  valorobjetivo <- list()
  status_vector <- list()
  j=j+1
  valorobjetivo_aux <- c()
  betas_aux <-c()
  status_vector_aux <- c()
  A_aux <- cbind(rep(0,ncol(Xtrain)-1), diag(ncol(Xtrain)-1)) #Tendra que cambiarse segun el tipo de problema
  A <- as.matrix(cbind(A_aux, -diag(nrow(A_aux)), diag(nrow(A_aux))))
  rhs        <- 0
  sense      <- '='
  objcon <- (1/nrow(Xtrain))*t(ytrain)%*%ytrain
  matrizcompuesta <- (1/nrow(Xtrain))*as.matrix(rbind(cbind(t(Xtrain)%*%Xtrain+lambda*(1-alpha)*(t(A_aux)%*%A_aux),matrix(rep(0,ncol(Xtrain)*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(Xtrain)+2*nrow(A))),nrow=2*nrow(A))))
  obj <- c(rep(0, ncol(Xtrain)), rep(alpha*lambda, 2*nrow(A))) - (2/nrow(Xtrain))*t(ytrain)%*%cbind(Xtrain,matrix(rep(0, 2*nrow(A)*nrow(Xtrain)),ncol=2*nrow(A)))
  Q <- matrizcompuesta
  vtype <- 'C'
  modelsense <- 'min'
  modelname <- 'CLassoMatrizId'
  lb <- c(rep(-Inf,ncol(Xtrain)), rep(0,2*nrow(A)))
  quadcon <- list()
  Cuadratica <- list()
  constraints_datasets <- list(list())
  constraints_datasets[[1]] <- cbind(y_1train,X_1train)
  constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
  constraints_datasets[[2]] <- cbind(y_2train,X_2train)
  constraints_datasets[[2]] <- as.matrix(constraints_datasets[[2]])
  constraints_datasets[[3]] <- cbind(y_3train,X_3train)
  constraints_datasets[[3]] <- as.matrix(constraints_datasets[[3]])
  constraint_number <- length(constraints_datasets)
  
  for (i in 1:constraint_number){
    matrizcompuestarestriccion <- as.matrix(rbind(cbind(t(constraints_datasets[[i]][,-1])%*%constraints_datasets[[i]][,-1],matrix(rep(0,ncol(constraints_datasets[[i]][,-1])*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(constraints_datasets[[i]][,-1])+2*nrow(A))),nrow=2*nrow(A))))
    Cuadratica$Qc[[i]]= (1/nrow(constraints_datasets[[i]][,-1]))*matrizcompuestarestriccion
    Cuadratica$q[[i]] = -(2/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%cbind(constraints_datasets[[i]][,-1],matrix(rep(0, 2*nrow(A)*nrow(constraints_datasets[[i]][,-1])),ncol=2*nrow(A)))
    Cuadratica$rhs[[i]] =  f[i] - (1/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%constraints_datasets[[i]][,1]
    Cuadratica$sense[[i]] <- '<'
    quadcon[[i]] <- list(Qc =  Cuadratica$Qc[[i]], q =  as.numeric(Cuadratica$q[[i]]), rhs= Cuadratica$rhs[[i]], sense = Cuadratica$sense[[i]])}
  
  model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
  result <- gurobi(model)
  status_vector_aux <- c(status_vector_aux, result$status)
  valorobjetivo_aux1 <- result$objval
  betas_aux <- result$x[1:ncol(Xtrain)]
  Xtest <- samples[folds[[g]],-c(22,23)]
  ytest <-samples[folds[[g]],22]
  residuos <- as.matrix(Xtest)%*%betas_aux - ytest
  MSE_aux <- sum(residuos^2)/length(residuos)
  MSE <- c(MSE,MSE_aux)
  betas[,which(gammas == gamma)] <- betas_aux
}



i=4

gammas[i]
length(which(betas[,i]>1e-7))
residuos <- as.matrix(Xtrain)%*%betas[,i] - ytrain
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_1)%*%betas[,i] - y_1
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_2)%*%betas[,i] - y_2
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_3)%*%betas[,i] - y_3
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_4)%*%betas[,i] - y_4
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_5)%*%betas[,i] - y_5
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_6)%*%betas[,i] - y_6
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_7)%*%betas[,i] - y_7
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux
residuos <- as.matrix(X_8)%*%betas[,i] - y_8
MSE_aux <- sum(residuos^2)/length(residuos)
MSE_aux



# Instalar glmnet si no está instalado
# install.packages("glmnet")

library(glmnet)

# Definir un vector de lambdas que quieres evaluar
lambda_seq <- vector_lambdas/2

# Entrenar modelo Lasso con validación cruzada
set.seed(123)  # para reproducibilidad
cv_fit <- cv.glmnet(X[,-1], y, alpha = 0.25, lambda = lambda_seq, standardize = TRUE)

# Obtener el mejor lambda que minimiza el MSE
best_lambda <- cv_fit$lambda.min

# Ajustar el modelo Lasso con el mejor lambda
final_model <- glmnet(X[,-1], y, alpha = 0.25, lambda = best_lambda, standardize = TRUE)
final_model$lambda
# Obtener el vector de coeficientes incluyendo el intercepto
beta <- as.vector(coef(final_model))

# Mostrar el vector beta
print(beta)


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
