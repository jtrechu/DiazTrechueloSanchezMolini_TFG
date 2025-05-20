# Function to solve single-constraint version of tau minimization problem
# Arguments:
#   constraints_datasets: List containing one constraint dataset
#   ECM_star: Single target error value
tau_min_function_2 <-  function(constraints_datasets, ECM_star){
  # Convert first dataset to matrix format
  constraints_datasets[[1]] <- as.matrix(constraints_datasets[[1]])
  
  # Initialize empty linear constraint (1 row of zeros)
  A = matrix(rep(0,ncol(constraints_datasets[[1]][,-1])+1), nrow = 1)
  rhs        <- 0  # Right-hand side
  sense      <- '='  # Equality constraint
  
  # Initialize quadratic objective matrix (all zeros)
  Q <- matrix(0, ncol= ncol(constraints_datasets[[1]][,-1])+1, nrow=ncol(constraints_datasets[[1]][,-1])+1)
  
  # Linear objective (minimize last variable)
  obj <- c(rep(0,ncol(constraints_datasets[[1]][,-1])),1)
  objcon <- 0  # Constant term in objective
  vtype <- 'C'  # Continuous variables
  modelsense <- 'min'  # Minimization problem
  modelname <- 'tau_min_problem'  # Model identifier
  
  # Variable bounds (unbounded except last >= 0)
  lb <- c(rep(-Inf, ncol(constraints_datasets[[1]][,-1])), 0)
  
  # Setup quadratic constraint components
  quadcon <- list()  # Will store final constraint
  Cuadratica <- list()  # Temporary storage
  
  # Construct quadratic constraint matrix X'X where X is data without first column
  matrizcompuestarestriccion <- as.matrix(t(cbind(constraints_datasets[[1]][,-1],rep(0, nrow(constraints_datasets[[1]][,-1]))))%*%cbind(constraints_datasets[[1]][,-1],rep(0, nrow(constraints_datasets[[1]][,-1]))))
  
  # Store scaled quadratic component (1/n * X'X)
  Cuadratica$Qc[[1]]= (1/nrow(constraints_datasets[[1]][,-1]))*matrizcompuestarestriccion
  
  # Construct linear component: -2/n * y'X - [0,...,0, ECM_star]
  q_temp <- -(2/nrow(constraints_datasets[[1]][,-1]))*t(constraints_datasets[[1]][,1])%*%cbind(constraints_datasets[[1]][,-1],rep(0, nrow(constraints_datasets[[1]][,-1]))) - c(rep(0,ncol(constraints_datasets[[1]][,-1])), ECM_star[1])
  Cuadratica$q[[1]] = as.numeric(q_temp)
  
  # NOTE: There appears to be a potential typo using 'i' instead of '1' in ECM_star[i]
  Cuadratica$rhs[[1]] = ECM_star[i] - (1/nrow(constraints_datasets[[1]][,-1]))*t(constraints_datasets[[1]][,1])%*%constraints_datasets[[1]][,1]    
  Cuadratica$sense[[1]] <- '<'  # Inequality constraint
  
  # Store complete quadratic constraint
  quadcon[[1]] <- list(Qc =  Cuadratica$Qc[[1]], q =  Cuadratica$q[[1]], rhs= Cuadratica$rhs[[1]], sense = Cuadratica$sense[[1]])
  
  # Build complete optimization model
  model <- list(A=A,rhs=rhs,sense=sense,obj=obj, Q=Q, vtype=vtype,modelsense=modelsense, lb=lb, objcon=objcon, quadcon=quadcon)
  
  # Solve using Gurobi optimizer
  result <- gurobi(model)
  
  # Return optimal objective value (minimum tau)
  tau_min_function <- result$objval
}