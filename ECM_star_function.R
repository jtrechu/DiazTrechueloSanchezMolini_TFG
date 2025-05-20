
ECM_star_function <- function(X,y){
  
  A = matrix(rep(0,ncol(X)), nrow = 1)
  rhs        <- 0
  sense      <- '='
  objcon <- (1/nrow(X))*t(y)%*%y
  matrizcompuesta <- (1/nrow(X))*as.matrix(t(X)%*%X)
  Q <- matrizcompuesta
  obj <- -(2/nrow(X))*t(y)%*%X
  vtype <- 'C'
  modelsense <- 'min'
  modelname <- 'ECM_star_problem'
  lb <- c(rep(-Inf,ncol(X)))
  
  model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense,objcon=objcon, lb=lb)
  
  result <- gurobi(model)
  
  ECM_star_function <- c(result$status, result$objval, result$x)
  
}