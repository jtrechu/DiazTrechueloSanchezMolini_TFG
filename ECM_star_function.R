ECM_star_function <- function(X, y) {
  # Initialize the linear constraint matrix A with zeros
  # Here, A is a 1-row matrix with number of columns equal to number of predictors
  A <- matrix(rep(0, ncol(X)), nrow = 1)
  
  # Right-hand side of the constraint (equals zero)
  rhs <- 0
  
  # Sense of the constraint: equality '='
  sense <- '='
  
  # Constant term in the objective function: (1/n) * y'y
  # This is the constant part of the quadratic objective function (doesn't depend on beta)
  objcon <- (1 / nrow(X)) * t(y) %*% y
  
  # Quadratic coefficient matrix for the objective function: (1/n) * X'X
  matrizcompuesta <- (1 / nrow(X)) * as.matrix(t(X) %*% X)
  Q <- matrizcompuesta
  
  # Linear coefficient vector for the objective function: -(2/n) * y'X
  # This term comes from expanding the squared error loss
  obj <- -(2 / nrow(X)) * t(y) %*% X
  
  # Variable type: continuous ('C')
  vtype <- 'C'
  
  # Model sense: minimize the objective function
  modelsense <- 'min'
  
  # Model name (optional)
  modelname <- 'ECM_star_problem'
  
  # Lower bounds on variables (no lower bound, i.e. -Inf)
  lb <- rep(-Inf, ncol(X))
  
  # Construct the model list for Gurobi solver
  model <- list(
    A = A,
    rhs = rhs,
    sense = sense,
    obj = obj,
    Q = Q,
    vtype = vtype,
    modelsense = modelsense,
    objcon = objcon,
    lb = lb
  )
  
  # Solve the quadratic programming problem with Gurobi
  result <- gurobi(model)
  
  # Return a vector with the optimization status, objective value, and solution vector x (beta)
  ECM_star_function <- c(result$status, result$objval, result$x)
}
