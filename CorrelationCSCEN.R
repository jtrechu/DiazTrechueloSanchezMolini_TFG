
options(digits = 11)
library(gurobi)
source('ECM_star_function.R')
source('tau_min_function_2.R')
source('tau_max_function.R')
source('CorrelationDataSet.R')
alpha<-0.25
p<-.99
vector_alphas <- c(alpha)
ext_lambdas <- 30
grosor_rejilla_lambdas <- 0.1
vector_lambdas <- seq(0,ext_lambdas,by=grosor_rejilla_lambdas)
grosor_rejilla_taus <-0.1
df <- generate_dataset(p,100,80)
df <- cbind(Intercept=1,df)
X <- df[,-c(2,13)]
X <- as.matrix(X)
y <-  df[,2]
X_2 <- as.matrix(df[which(df[,13]==2),-c(2,13)])
y_2 <- df[which(df[,13]==2),2]
##Computation of MSE*_l:
constraints_datasets <- list(list())
constraints_datasets[[1]] <- cbind(y_2,X_2)
constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
constraint_number <- length(constraints_datasets)
Sol_OLS_constraints <- c()
f0_aux <- ECM_star_function(as.matrix(constraints_datasets[[1]][,-1]), constraints_datasets[[1]][,1])
ECM_star <- as.numeric(f0_aux[2])
Sol_OLS_constraints <-as.numeric(f0_aux[3:(ncol(constraints_datasets[[1]][,-1])+2)])
vector_tauss<-c(50)
if (is.null(vector_tauss)) {
  #Computation of tau_min (no dependence on hyperparameters)
  tau_min <- tau_min_function_2(constraints_datasets, ECM_star) ##REVISAR
  tau_max_lambda_vector <- tau_max(constraints_datasets,vector_lambdas,alpha)
  vector_tauss <-  seq(tau_min, max(tau_max_lambda_vector)+2, by=grosor_rejilla_taus)}
vector_tauss
j=0
valorobjetivo <- list()
betas <- list()
status_vector <- list()

for (alpha in vector_alphas){
  for (lambda in vector_lambdas){
    j=j+1
    valorobjetivo_aux <- c()
    betas_aux <-c()
    status_vector_aux <- c()
    vector_tauss <- c(50)
    for (tau in vector_tauss){
      
      ## Funcion Objetivo
      A_aux <- cbind(rep(0,ncol(X)-1), diag(ncol(X)-1)) #Tendra que cambiarse segun el tipo de problema
      A <- as.matrix(cbind(A_aux, -diag(nrow(A_aux)), diag(nrow(A_aux))))
      rhs        <- 0
      sense      <- '='
      objcon <- (1/nrow(X))*t(y)%*%y
      matrizcompuesta <- (1/nrow(X))*as.matrix(rbind(cbind(t(X)%*%X+lambda*(1-alpha)*(t(A_aux)%*%A_aux),matrix(rep(0,ncol(X)*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(X)+2*nrow(A))),nrow=2*nrow(A))))
      obj <- c(rep(0, ncol(X)), rep(alpha*lambda, 2*nrow(A))) - (2/nrow(X))*t(y)%*%cbind(X,matrix(rep(0, 2*nrow(A)*nrow(X)),ncol=2*nrow(A)))
      Q <- matrizcompuesta
      vtype <- 'C'
      modelsense <- 'min'
      modelname <- 'CLassoMatrizId'
      lb <- c(rep(-Inf,ncol(X)), rep(0,2*nrow(A)))
      quadcon <- list()
      Cuadratica <- list()
      for (i in 1:length(constraints_datasets)){
        matrizcompuestarestriccion <- as.matrix(rbind(cbind(t(constraints_datasets[[i]][,-1])%*%constraints_datasets[[i]][,-1],matrix(rep(0,ncol(constraints_datasets[[i]][,-1])*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(constraints_datasets[[i]][,-1])+2*nrow(A))),nrow=2*nrow(A))))
        Cuadratica$Qc[[i]]= (1/nrow(constraints_datasets[[i]][,-1]))*matrizcompuestarestriccion
        Cuadratica$q[[i]] = -(2/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%cbind(constraints_datasets[[i]][,-1],matrix(rep(0, 2*nrow(A)*nrow(constraints_datasets[[i]][,-1])),ncol=2*nrow(A)))
        Cuadratica$rhs[[i]] =  (1+tau)*ECM_star[i] - (1/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%constraints_datasets[[i]][,1]
        Cuadratica$sense[[i]] <- '<'
        
        quadcon[[i]] <- list(Qc =  Cuadratica$Qc[[i]], q =  as.numeric(Cuadratica$q[[i]]), rhs= Cuadratica$rhs[[i]], sense = Cuadratica$sense[[i]])
      }
      
      model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
      
      result <- gurobi(model)
      status_vector_aux <- c(status_vector_aux, result$status)
      
      
      valorobjetivo_aux1 <- result$objval
      betas_aux1 <- result$x[1:ncol(X)]
      
      valorobjetivo_aux <- c(valorobjetivo_aux, valorobjetivo_aux1)
      betas_aux <- rbind(betas_aux, betas_aux1)
    }
    
    valorobjetivo[[j]] <- valorobjetivo_aux
    #betas[[j]] <- rbind(matrix(logical(0), nrow=(length(vector_tauss)- length(tauss)), ncol=ncol(X)),betas_aux)
    #betas[[j]] <- rbind(betas_aux, matrix(logical(0), nrow=(length(vector_tauss)- length(tauss)), ncol=ncol(X)))
    betas[[j]] <- betas_aux
    status_vector[[j]] <- status_vector_aux
  }}

lambda <- seq(0,30,by=0.100)
beta_matrix <- do.call(rbind, betas)
beta_matrix <- cbind(lambda,beta_matrix[,c(2,3,4)])
colnames(beta_matrix) <- c("lambda",paste0("beta", 1:3))
df <- as.data.frame(beta_matrix)

DifEN3 <- abs(df$beta1-df$beta2)
##ALPHA=0.5
options(digits = 11)
library(gurobi)
path_directory = "C:/Users/jtrec/OneDrive/Desktop/Codigo Jaime/"
source(paste(path_directory, 'ECM_star_function.R', sep=""))
source(paste(path_directory, 'tau_min_function_2.R', sep=""))
source(paste(path_directory, 'tau_max_function.R', sep=""))
source(paste(path_directory, 'CorrelationDataSet.R', sep=""))
alpha<-0.5
p<-.99
vector_alphas <- c(alpha)
ext_lambdas <- 30
grosor_rejilla_lambdas <- 0.1
vector_lambdas <- seq(0,ext_lambdas,by=grosor_rejilla_lambdas)
grosor_rejilla_taus <-0.1
df <- generate_dataset(p,100,80)
df <- cbind(Intercept=1,df)
X <- df[,-c(2,13)]
X <- as.matrix(X)
y <-  df[,2]
X_2 <- as.matrix(df[which(df[,13]==2),-c(2,13)])
y_2 <- df[which(df[,13]==2),2]
##Computation of MSE*_l:
constraints_datasets <- list(list())
constraints_datasets[[1]] <- cbind(y_2,X_2)
constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
constraint_number <- length(constraints_datasets)
Sol_OLS_constraints <- c()
f0_aux <- ECM_star_function(as.matrix(constraints_datasets[[1]][,-1]), constraints_datasets[[1]][,1])
ECM_star <- as.numeric(f0_aux[2])
Sol_OLS_constraints <-as.numeric(f0_aux[3:(ncol(constraints_datasets[[1]][,-1])+2)])
vector_tauss<-c(50)
if (is.null(vector_tauss)) {
  #Computation of tau_min (no dependence on hyperparameters)
  tau_min <- tau_min_function_2(constraints_datasets, ECM_star) ##REVISAR
  tau_max_lambda_vector <- tau_max(constraints_datasets,vector_lambdas,alpha)
  vector_tauss <-  seq(tau_min, max(tau_max_lambda_vector)+2, by=grosor_rejilla_taus)}
vector_tauss
j=0
valorobjetivo <- list()
betas <- list()
status_vector <- list()

for (alpha in vector_alphas){
  for (lambda in vector_lambdas){
    j=j+1
    valorobjetivo_aux <- c()
    betas_aux <-c()
    status_vector_aux <- c()
    vector_tauss <- c(50)
    for (tau in vector_tauss){
      
      ## Funcion Objetivo
      A_aux <- cbind(rep(0,ncol(X)-1), diag(ncol(X)-1)) #Tendra que cambiarse segun el tipo de problema
      A <- as.matrix(cbind(A_aux, -diag(nrow(A_aux)), diag(nrow(A_aux))))
      rhs        <- 0
      sense      <- '='
      objcon <- (1/nrow(X))*t(y)%*%y
      matrizcompuesta <- (1/nrow(X))*as.matrix(rbind(cbind(t(X)%*%X+lambda*(1-alpha)*(t(A_aux)%*%A_aux),matrix(rep(0,ncol(X)*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(X)+2*nrow(A))),nrow=2*nrow(A))))
      obj <- c(rep(0, ncol(X)), rep(alpha*lambda, 2*nrow(A))) - (2/nrow(X))*t(y)%*%cbind(X,matrix(rep(0, 2*nrow(A)*nrow(X)),ncol=2*nrow(A)))
      Q <- matrizcompuesta
      vtype <- 'C'
      modelsense <- 'min'
      modelname <- 'CLassoMatrizId'
      lb <- c(rep(-Inf,ncol(X)), rep(0,2*nrow(A)))
      quadcon <- list()
      Cuadratica <- list()
      for (i in 1:length(constraints_datasets)){
        matrizcompuestarestriccion <- as.matrix(rbind(cbind(t(constraints_datasets[[i]][,-1])%*%constraints_datasets[[i]][,-1],matrix(rep(0,ncol(constraints_datasets[[i]][,-1])*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(constraints_datasets[[i]][,-1])+2*nrow(A))),nrow=2*nrow(A))))
        Cuadratica$Qc[[i]]= (1/nrow(constraints_datasets[[i]][,-1]))*matrizcompuestarestriccion
        Cuadratica$q[[i]] = -(2/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%cbind(constraints_datasets[[i]][,-1],matrix(rep(0, 2*nrow(A)*nrow(constraints_datasets[[i]][,-1])),ncol=2*nrow(A)))
        Cuadratica$rhs[[i]] =  (1+tau)*ECM_star[i] - (1/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%constraints_datasets[[i]][,1]
        Cuadratica$sense[[i]] <- '<'
        
        quadcon[[i]] <- list(Qc =  Cuadratica$Qc[[i]], q =  as.numeric(Cuadratica$q[[i]]), rhs= Cuadratica$rhs[[i]], sense = Cuadratica$sense[[i]])
      }
      
      model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
      
      result <- gurobi(model)
      status_vector_aux <- c(status_vector_aux, result$status)
      
      
      valorobjetivo_aux1 <- result$objval
      betas_aux1 <- result$x[1:ncol(X)]
      
      valorobjetivo_aux <- c(valorobjetivo_aux, valorobjetivo_aux1)
      betas_aux <- rbind(betas_aux, betas_aux1)
    }
    
    valorobjetivo[[j]] <- valorobjetivo_aux
    #betas[[j]] <- rbind(matrix(logical(0), nrow=(length(vector_tauss)- length(tauss)), ncol=ncol(X)),betas_aux)
    #betas[[j]] <- rbind(betas_aux, matrix(logical(0), nrow=(length(vector_tauss)- length(tauss)), ncol=ncol(X)))
    betas[[j]] <- betas_aux
    status_vector[[j]] <- status_vector_aux
  }}


library(ggplot2)
library(latex2exp)
library(tidyr)
library(dplyr)

lambda <- seq(0,30,by=0.100)
beta_matrix <- do.call(rbind, betas)
beta_matrix <- cbind(lambda,beta_matrix[,c(2,3,4)])
colnames(beta_matrix) <- c("lambda",paste0("beta", 1:3))
df <- as.data.frame(beta_matrix)

DifEN <- abs(df$beta1-df$beta2)

##ALPHA=0.75
options(digits = 11)
library(gurobi)
path_directory = "C:/Users/jtrec/OneDrive/Desktop/Codigo Jaime/"
source(paste(path_directory, 'ECM_star_function.R', sep=""))
source(paste(path_directory, 'tau_min_function_2.R', sep=""))
source(paste(path_directory, 'tau_max_function.R', sep=""))
source(paste(path_directory, 'CorrelationDataSet.R', sep=""))
alpha<-0.75
p<-.99
vector_alphas <- c(alpha)
ext_lambdas <- 30
grosor_rejilla_lambdas <- 0.1
vector_lambdas <- seq(0,ext_lambdas,by=grosor_rejilla_lambdas)
grosor_rejilla_taus <-0.1
df <- generate_dataset(p,100,80)
df <- cbind(Intercept=1,df)
X <- df[,-c(2,13)]
X <- as.matrix(X)
y <-  df[,2]
X_2 <- as.matrix(df[which(df[,13]==2),-c(2,13)])
y_2 <- df[which(df[,13]==2),2]
##Computation of MSE*_l:
constraints_datasets <- list(list())
constraints_datasets[[1]] <- cbind(y_2,X_2)
constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
constraint_number <- length(constraints_datasets)
Sol_OLS_constraints <- c()
f0_aux <- ECM_star_function(as.matrix(constraints_datasets[[1]][,-1]), constraints_datasets[[1]][,1])
ECM_star <- as.numeric(f0_aux[2])
Sol_OLS_constraints <-as.numeric(f0_aux[3:(ncol(constraints_datasets[[1]][,-1])+2)])
vector_tauss<-c(50)
if (is.null(vector_tauss)) {
  #Computation of tau_min (no dependence on hyperparameters)
  tau_min <- tau_min_function_2(constraints_datasets, ECM_star) ##REVISAR
  tau_max_lambda_vector <- tau_max(constraints_datasets,vector_lambdas,alpha)
  vector_tauss <-  seq(tau_min, max(tau_max_lambda_vector)+2, by=grosor_rejilla_taus)}
vector_tauss
j=0
valorobjetivo <- list()
betas <- list()
status_vector <- list()

for (alpha in vector_alphas){
  for (lambda in vector_lambdas){
    j=j+1
    valorobjetivo_aux <- c()
    betas_aux <-c()
    status_vector_aux <- c()
    vector_tauss <- c(50)
    for (tau in vector_tauss){
      
      ## Funcion Objetivo
      A_aux <- cbind(rep(0,ncol(X)-1), diag(ncol(X)-1)) #Tendra que cambiarse segun el tipo de problema
      A <- as.matrix(cbind(A_aux, -diag(nrow(A_aux)), diag(nrow(A_aux))))
      rhs        <- 0
      sense      <- '='
      objcon <- (1/nrow(X))*t(y)%*%y
      matrizcompuesta <- (1/nrow(X))*as.matrix(rbind(cbind(t(X)%*%X+lambda*(1-alpha)*(t(A_aux)%*%A_aux),matrix(rep(0,ncol(X)*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(X)+2*nrow(A))),nrow=2*nrow(A))))
      obj <- c(rep(0, ncol(X)), rep(alpha*lambda, 2*nrow(A))) - (2/nrow(X))*t(y)%*%cbind(X,matrix(rep(0, 2*nrow(A)*nrow(X)),ncol=2*nrow(A)))
      Q <- matrizcompuesta
      vtype <- 'C'
      modelsense <- 'min'
      modelname <- 'CLassoMatrizId'
      lb <- c(rep(-Inf,ncol(X)), rep(0,2*nrow(A)))
      quadcon <- list()
      Cuadratica <- list()
      for (i in 1:length(constraints_datasets)){
        matrizcompuestarestriccion <- as.matrix(rbind(cbind(t(constraints_datasets[[i]][,-1])%*%constraints_datasets[[i]][,-1],matrix(rep(0,ncol(constraints_datasets[[i]][,-1])*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(constraints_datasets[[i]][,-1])+2*nrow(A))),nrow=2*nrow(A))))
        Cuadratica$Qc[[i]]= (1/nrow(constraints_datasets[[i]][,-1]))*matrizcompuestarestriccion
        Cuadratica$q[[i]] = -(2/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%cbind(constraints_datasets[[i]][,-1],matrix(rep(0, 2*nrow(A)*nrow(constraints_datasets[[i]][,-1])),ncol=2*nrow(A)))
        Cuadratica$rhs[[i]] =  (1+tau)*ECM_star[i] - (1/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%constraints_datasets[[i]][,1]
        Cuadratica$sense[[i]] <- '<'
        
        quadcon[[i]] <- list(Qc =  Cuadratica$Qc[[i]], q =  as.numeric(Cuadratica$q[[i]]), rhs= Cuadratica$rhs[[i]], sense = Cuadratica$sense[[i]])
      }
      
      model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
      
      result <- gurobi(model)
      status_vector_aux <- c(status_vector_aux, result$status)
      
      
      valorobjetivo_aux1 <- result$objval
      betas_aux1 <- result$x[1:ncol(X)]
      
      valorobjetivo_aux <- c(valorobjetivo_aux, valorobjetivo_aux1)
      betas_aux <- rbind(betas_aux, betas_aux1)
    }
    
    valorobjetivo[[j]] <- valorobjetivo_aux
    #betas[[j]] <- rbind(matrix(logical(0), nrow=(length(vector_tauss)- length(tauss)), ncol=ncol(X)),betas_aux)
    #betas[[j]] <- rbind(betas_aux, matrix(logical(0), nrow=(length(vector_tauss)- length(tauss)), ncol=ncol(X)))
    betas[[j]] <- betas_aux
    status_vector[[j]] <- status_vector_aux
  }}
lambda <- seq(0,30,by=0.100)
beta_matrix <- do.call(rbind, betas)
beta_matrix <- cbind(lambda,beta_matrix[,c(2,3,4)])
colnames(beta_matrix) <- c("lambda",paste0("beta", 1:3))
df <- as.data.frame(beta_matrix)

DifEN2 <- abs(df$beta1-df$beta2)

##ALPHA=1
options(digits = 11)
library(gurobi)
path_directory = "C:/Users/jtrec/OneDrive/Desktop/Codigo Jaime/"
source(paste(path_directory, 'ECM_star_function.R', sep=""))
source(paste(path_directory, 'tau_min_function_2.R', sep=""))
source(paste(path_directory, 'tau_max_function.R', sep=""))
source(paste(path_directory, 'CorrelationDataSet.R', sep=""))
alpha<-1
p<-.99
vector_alphas <- c(alpha)
ext_lambdas <- 30
grosor_rejilla_lambdas <- 0.1
vector_lambdas <- seq(0,ext_lambdas,by=grosor_rejilla_lambdas)
grosor_rejilla_taus <-0.1
df <- generate_dataset(p,100,80)
df <- cbind(Intercept=1,df)
X <- df[,-c(2,13)]
X <- as.matrix(X)
y <-  df[,2]
X_2 <- as.matrix(df[which(df[,13]==2),-c(2,13)])
y_2 <- df[which(df[,13]==2),2]
##Computation of MSE*_l:
constraints_datasets <- list(list())
constraints_datasets[[1]] <- cbind(y_2,X_2)
constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
constraint_number <- length(constraints_datasets)
Sol_OLS_constraints <- c()
f0_aux <- ECM_star_function(as.matrix(constraints_datasets[[1]][,-1]), constraints_datasets[[1]][,1])
ECM_star <- as.numeric(f0_aux[2])
Sol_OLS_constraints <-as.numeric(f0_aux[3:(ncol(constraints_datasets[[1]][,-1])+2)])
vector_tauss<-c(50)
if (is.null(vector_tauss)) {
  #Computation of tau_min (no dependence on hyperparameters)
  tau_min <- tau_min_function_2(constraints_datasets, ECM_star) ##REVISAR
  tau_max_lambda_vector <- tau_max(constraints_datasets,vector_lambdas,alpha)
  vector_tauss <-  seq(tau_min, max(tau_max_lambda_vector)+2, by=grosor_rejilla_taus)}
vector_tauss
j=0
valorobjetivo <- list()
betas <- list()
status_vector <- list()

for (alpha in vector_alphas){
  for (lambda in vector_lambdas){
    j=j+1
    valorobjetivo_aux <- c()
    betas_aux <-c()
    status_vector_aux <- c()
    vector_tauss <- c(50)
    for (tau in vector_tauss){
      
      ## Funcion Objetivo
      A_aux <- cbind(rep(0,ncol(X)-1), diag(ncol(X)-1)) #Tendra que cambiarse segun el tipo de problema
      A <- as.matrix(cbind(A_aux, -diag(nrow(A_aux)), diag(nrow(A_aux))))
      rhs        <- 0
      sense      <- '='
      objcon <- (1/nrow(X))*t(y)%*%y
      matrizcompuesta <- (1/nrow(X))*as.matrix(rbind(cbind(t(X)%*%X+lambda*(1-alpha)*(t(A_aux)%*%A_aux),matrix(rep(0,ncol(X)*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(X)+2*nrow(A))),nrow=2*nrow(A))))
      obj <- c(rep(0, ncol(X)), rep(alpha*lambda, 2*nrow(A))) - (2/nrow(X))*t(y)%*%cbind(X,matrix(rep(0, 2*nrow(A)*nrow(X)),ncol=2*nrow(A)))
      Q <- matrizcompuesta
      vtype <- 'C'
      modelsense <- 'min'
      modelname <- 'CLassoMatrizId'
      lb <- c(rep(-Inf,ncol(X)), rep(0,2*nrow(A)))
      quadcon <- list()
      Cuadratica <- list()
      for (i in 1:length(constraints_datasets)){
        matrizcompuestarestriccion <- as.matrix(rbind(cbind(t(constraints_datasets[[i]][,-1])%*%constraints_datasets[[i]][,-1],matrix(rep(0,ncol(constraints_datasets[[i]][,-1])*2*nrow(A)), ncol=2*nrow(A))),matrix(rep(0, 2*nrow(A)*(ncol(constraints_datasets[[i]][,-1])+2*nrow(A))),nrow=2*nrow(A))))
        Cuadratica$Qc[[i]]= (1/nrow(constraints_datasets[[i]][,-1]))*matrizcompuestarestriccion
        Cuadratica$q[[i]] = -(2/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%cbind(constraints_datasets[[i]][,-1],matrix(rep(0, 2*nrow(A)*nrow(constraints_datasets[[i]][,-1])),ncol=2*nrow(A)))
        Cuadratica$rhs[[i]] =  (1+tau)*ECM_star[i] - (1/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%constraints_datasets[[i]][,1]
        Cuadratica$sense[[i]] <- '<'
        
        quadcon[[i]] <- list(Qc =  Cuadratica$Qc[[i]], q =  as.numeric(Cuadratica$q[[i]]), rhs= Cuadratica$rhs[[i]], sense = Cuadratica$sense[[i]])
      }
      
      model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
      
      result <- gurobi(model)
      status_vector_aux <- c(status_vector_aux, result$status)
      
      
      valorobjetivo_aux1 <- result$objval
      betas_aux1 <- result$x[1:ncol(X)]
      
      valorobjetivo_aux <- c(valorobjetivo_aux, valorobjetivo_aux1)
      betas_aux <- rbind(betas_aux, betas_aux1)
    }
    
    valorobjetivo[[j]] <- valorobjetivo_aux
    #betas[[j]] <- rbind(matrix(logical(0), nrow=(length(vector_tauss)- length(tauss)), ncol=ncol(X)),betas_aux)
    #betas[[j]] <- rbind(betas_aux, matrix(logical(0), nrow=(length(vector_tauss)- length(tauss)), ncol=ncol(X)))
    betas[[j]] <- betas_aux
    status_vector[[j]] <- status_vector_aux
  }}


library(ggplot2)
library(latex2exp)
library(tidyr)
library(dplyr)

lambda <- seq(0,30,by=0.100)
beta_matrix <- do.call(rbind, betas)
beta_matrix <- cbind(lambda,beta_matrix[,c(2,3,4)])
colnames(beta_matrix) <- c("lambda",paste0("beta", 1:3))
df <- as.data.frame(beta_matrix)

DifLasso <- abs(df$beta1-df$beta2)

library(ggplot2)
library(latex2exp)
library(tidyr)
library(dplyr)

lambda <- seq(0,30,by=0.100)
beta_matrix <- do.call(rbind, betas)
beta_matrix <- cbind(lambda,beta_matrix[,c(2,3,4)])
colnames(beta_matrix) <- c("lambda",paste0("beta", 1:3))
df <- as.data.frame(beta_matrix)


library(ggplot2)
library(latex2exp)
library(tidyr)
library(dplyr)

lambda <- seq(0,30,by=0.100)
beta_matrix <- do.call(rbind, betas)
beta_matrix <- cbind(lambda,beta_matrix[,c(2,3,4)])
colnames(beta_matrix) <- c("lambda",paste0("beta", 1:3))
df <- as.data.frame(beta_matrix)



diferencias <- data.frame(lambda=lambda, Lasso = DifLasso, EN05 = DifEN, EN075 = DifEN2, EN025 = DifEN3)

diferencias_long <- pivot_longer(diferencias, cols = c("Lasso", "EN025","EN05", "EN075"),
                                 names_to = "Method", values_to = "Difference")

# Plot
ggplot(diferencias_long, aes(x = lambda, y = Difference, color = Method)) +
  geom_line(size = 1) +
  labs(title = "Comparison of CSCLasso and CSCEN",
       x = expression(lambda),
       y = TeX("$|\\beta_1-\\beta_2|$"),
       color = "Method") +
  scale_color_manual(values = c("Lasso" = "blue", "EN025" = "red", "EN05" = "orange", "EN075"="purple"),
                     labels = c("Lasso" = TeX("CSCLasso ($\\alpha=1$)"), "EN025" = TeX("CSCEN ($\\alpha=0.25$)"), "EN05" = TeX("CSCEN ($\\alpha=0.5$)"), "EN075" = TeX("CSCEN ($\\alpha=0.75$)")))+
  theme_minimal()




