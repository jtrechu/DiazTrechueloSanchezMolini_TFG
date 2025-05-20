tau_min_function <-  function(constraints_datasets, ECM_star){
  constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
  constraints_datasets[[2]] <- as.matrix(constraints_datasets[[2]])
  A = matrix(rep(0,ncol(constraints_datasets[[1]][,-1])+1), nrow = 1)
  rhs        <- 0
  sense      <- '='
  Q <- matrix(0, ncol= ncol(constraints_datasets[[1]][,-1])+1, nrow=ncol(constraints_datasets[[1]][,-1])+1)
  obj <- c(rep(0,ncol(constraints_datasets[[1]][,-1])),1)
  objcon <- 0
  vtype <- 'C'
  modelsense <- 'min'
  modelname <- 'tau_min_problem'
  lb <- c(rep(-Inf, ncol(constraints_datasets[[1]][,-1])), 0)
  
  quadcon <- list()
  Cuadratica <- list()
  for (i in 1:length(constraints_datasets)){
    matrizcompuestarestriccion <- as.matrix(t(cbind(constraints_datasets[[i]][,-1],rep(0, nrow(constraints_datasets[[i]][,-1]))))%*%cbind(constraints_datasets[[i]][,-1],rep(0, nrow(constraints_datasets[[i]][,-1]))))
    Cuadratica$Qc[[i]]= (1/nrow(constraints_datasets[[i]][,-1]))*matrizcompuestarestriccion
    q_temp <- -(2/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%cbind(constraints_datasets[[i]][,-1],rep(0, nrow(constraints_datasets[[i]][,-1]))) - c(rep(0,ncol(constraints_datasets[[i]][,-1])), ECM_star[i])
    Cuadratica$q[[i]] = as.numeric(q_temp)
    Cuadratica$rhs[[i]] = ECM_star[i] - (1/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%constraints_datasets[[i]][,1]
    Cuadratica$sense[[i]] <- '<'
    
    quadcon[[i]] <- list(Qc =  Cuadratica$Qc[[i]], q =  Cuadratica$q[[i]], rhs= Cuadratica$rhs[[i]], sense = Cuadratica$sense[[i]])
  }
  
  model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
  
  result <- gurobi(model)
  
  tau_min_function <- result$objval
}