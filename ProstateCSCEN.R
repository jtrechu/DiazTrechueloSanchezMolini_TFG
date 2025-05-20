CSCEN<- function(lambda,alpha,grosor_lambda=0.5,vector_tauss=NULL) {
  options(digits = 11)
  library(gurobi)
  source('ECM_star_function.R')
  source('tau_min_function.R')
  source('tau_max_function.R')
  
  vector_alphas <- c(alpha)
  ext_lambdas <- lambda
  grosor_rejilla_lambdas <- grosor_lambda
  vector_lambdas <- seq(0,ext_lambdas,by=grosor_rejilla_lambdas)
  grosor_rejilla_taus <-0.05
  library(gurobi)
  library(genridge)
  data(prostate)
  Prostate<-cbind(Intercept = 1, prostate[,-10])
  X <- Prostate[,-10]
  X <- as.matrix(X)
  y <-  Prostate[,10]
  X_menoresde65 <- Prostate[which(Prostate[,"age"]<65),-10]
  X_mayoresde65 <- Prostate[which(Prostate[,"age"]>=65),-10]
  y_menoresde65 <- Prostate[which(Prostate[,"age"]<65),10]
  y_mayoresde65 <- Prostate[which(Prostate[,"age"]>=65),10]
  
  ##Computation of MSE*_l:
  ECM_star <- c()
  constraints_datasets <- list(list())
  constraints_datasets[[1]] <- cbind(y_menoresde65,X_menoresde65)
  constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
  constraints_datasets[[2]] <- cbind(y_mayoresde65,X_mayoresde65)
  constraints_datasets[[2]] <- as.matrix(constraints_datasets[[2]])
  constraint_number <- length(constraints_datasets)
  Sol_OLS_constraints <- c()
  for (i in 1:constraint_number){
    f0_aux <- ECM_star_function(as.matrix(constraints_datasets[[i]][,-1]), constraints_datasets[[i]][,1])
    ECM_star_aux <- as.numeric(f0_aux[2])
    ECM_star <- c(ECM_star, ECM_star_aux)
    Sol_OLS_constraints <- rbind(Sol_OLS_constraints, as.numeric(f0_aux[3:(ncol(constraints_datasets[[i]][,-1])+2)]))
  }
  
  if (is.null(vector_tauss)) {
    #Computation of tau_min (no dependence on hyperparameters)
    tau_min <- tau_min_function(constraints_datasets, ECM_star)
    tau_max_lambda_vector <- tau_max(constraints_datasets,vector_lambdas,alpha)
    vector_tauss <-  seq(tau_min, max(tau_max_lambda_vector)+2, by=grosor_rejilla_taus)}
  
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
  return(list(betas=betas,taus=vector_tauss))
}


