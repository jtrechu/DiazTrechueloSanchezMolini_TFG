# Function to solve the optimization problem for minimum tau value
# Arguments:
#   constraints_datasets: List of constraint matrices
#   ECM_star: Vector of target error values
tau_min_function <-  function(constraints_datasets, ECM_star){
  # Convert datasets to matrices
  constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
  constraints_datasets[[2]] <- as.matrix(constraints_datasets[[2]])
  
  # Initialize linear constraint (empty)
  A = matrix(rep(0,ncol(constraints_datasets[[1]][,-1])+1), nrow = 1)
  rhs        <- 0
  sense      <- '='
  
  # Initialize quadratic objective (zero matrix)
  Q <- matrix(0, ncol= ncol(constraints_datasets[[1]][,-1])+1, nrow=ncol(constraints_datasets[[1]][,-1])+1)
  
  # Linear objective (minimize last variable)
  obj <- c(rep(0,ncol(constraints_datasets[[1]][,-1])),1)
  objcon <- 0
  vtype <- 'C'  # Continuous variables
  modelsense <- 'min'
  modelname <- 'tau_min_problem'
  
  # Variable bounds (last variable >= 0)
  lb <- c(rep(-Inf, ncol(constraints_datasets[[1]][,-1])), 0)
  
  # Prepare quadratic constraints
  quadcon <- list()
  Cuadratica <- list()
  
  for (i in 1:length(constraints_datasets)){
    # Construct quadratic constraint matrix
    matrizcompuestarestriccion <- as.matrix(t(cbind(constraints_datasets[[i]][,-1],rep(0, nrow(constraints_datasets[[i]][,-1]))))%*%cbind(constraints_datasets[[i]][,-1],rep(0, nrow(constraints_datasets[[i]][,-1]))))
    
    # Store scaled quadratic component
    Cuadratica$Qc[[i]]= (1/nrow(constraints_datasets[[i]][,-1]))*matrizcompuestarestriccion
    
    # Construct linear component of constraint
    q_temp <- -(2/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%cbind(constraints_datasets[[i]][,-1],rep(0, nrow(constraints_datasets[[i]][,-1]))) - c(rep(0,ncol(constraints_datasets[[i]][,-1])), ECM_star[i])
    Cuadratica$q[[i]] = as.numeric(q_temp)
    
    # Construct constraint RHS
    Cuadratica$rhs[[i]] = ECM_star[i] - (1/nrow(constraints_datasets[[i]][,-1]))*t(constraints_datasets[[i]][,1])%*%constraints_datasets[[i]][,1]
    Cuadratica$sense[[i]] <- '<'  # Inequality constraint
    
    # Store complete constraint
    quadcon[[i]] <- list(Qc =  Cuadratica$Qc[[i]], q =  Cuadratica$q[[i]], rhs= Cuadratica$rhs[[i]], sense = Cuadratica$sense[[i]])
  }
  
  # Build complete model
  model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
  
  # Solve with Gurobi
  result <- gurobi(model)
  
  # Return optimal value
  tau_min_function <- result$objval
}