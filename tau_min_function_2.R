tau_min_function_2 <-  function(constraints_datasets, ECM_star){
  constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
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
  matrizcompuestarestriccion <- as.matrix(t(cbind(constraints_datasets[[1]][,-1],rep(0, nrow(constraints_datasets[[1]][,-1]))))%*%cbind(constraints_datasets[[1]][,-1],rep(0, nrow(constraints_datasets[[1]][,-1]))))
  Cuadratica$Qc[[1]]= (1/nrow(constraints_datasets[[1]][,-1]))*matrizcompuestarestriccion
  q_temp <- -(2/nrow(constraints_datasets[[1]][,-1]))*t(constraints_datasets[[1]][,1])%*%cbind(constraints_datasets[[1]][,-1],rep(0, nrow(constraints_datasets[[1]][,-1]))) - c(rep(0,ncol(constraints_datasets[[1]][,-1])), ECM_star[1])
  Cuadratica$q[[1]] = as.numeric(q_temp)
  Cuadratica$rhs[[1]] = ECM_star[i] - (1/nrow(constraints_datasets[[1]][,-1]))*t(constraints_datasets[[1]][,1])%*%constraints_datasets[[1]][,1]    
  Cuadratica$sense[[1]] <- '<'
  quadcon[[1]] <- list(Qc =  Cuadratica$Qc[[1]], q =  Cuadratica$q[[1]], rhs= Cuadratica$rhs[[1]], sense = Cuadratica$sense[[1]])
  
  model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
  
  result <- gurobi(model)
  
  tau_min_function <- result$objval
}