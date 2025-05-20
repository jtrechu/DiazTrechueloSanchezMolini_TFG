CSCEN <- function(lambda, alpha, grosor_lambda = 0.5, vector_tauss = NULL) {
  # Set numeric printing precision
  options(digits = 11)
  
  # Load required libraries
  library(gurobi)
  library(genridge)
  
  # Source auxiliary functions for ECM star and tau computations
  source('ECM_star_function.R')
  source('tau_min_function.R')
  source('tau_max_function.R')
  
  # Initialize vectors of hyperparameters alpha and lambda
  vector_alphas <- c(alpha)
  ext_lambdas <- lambda
  
  # Define lambda grid step size (resolution)
  grosor_rejilla_lambdas <- grosor_lambda
  
  # Create sequence of lambda values from 0 up to ext_lambdas with step grosor_rejilla_lambdas
  vector_lambdas <- seq(0, ext_lambdas, by = grosor_rejilla_lambdas)
  
  # Define tau grid resolution
  grosor_rejilla_taus <- 0.05
  
  # Load prostate dataset and preprocess
  data(prostate)
  
  # Add intercept column, exclude response column (assumed to be column 10)
  Prostate <- cbind(Intercept = 1, prostate[, -10])
  X <- as.matrix(Prostate[, -10])   # Predictor matrix
  y <- Prostate[, 10]               # Response vector
  
  # Subset data by age groups for constraints later
  X_menoresde65 <- Prostate[which(Prostate[, "age"] < 65), -10]
  X_mayoresde65 <- Prostate[which(Prostate[, "age"] >= 65), -10]
  y_menoresde65 <- Prostate[which(Prostate[, "age"] < 65), 10]
  y_mayoresde65 <- Prostate[which(Prostate[, "age"] >= 65), 10]
  
  ## Compute ECM* (minimum squared errors) for the constraint datasets
  ECM_star <- c()
  constraints_datasets <- list(list())  # Initialize list to hold datasets for constraints
  
  # Add two constraint datasets combining y and X for each age group
  constraints_datasets[[1]] <- as.matrix(cbind(y_menoresde65, X_menoresde65))
  constraints_datasets[[2]] <- as.matrix(cbind(y_mayoresde65, X_mayoresde65))
  
  constraint_number <- length(constraints_datasets)
  Sol_OLS_constraints <- c()
  
  # Loop over constraint datasets to compute ECM_star and OLS solutions
  for (i in 1:constraint_number) {
    # Call ECM_star_function to get optimization status, objective, and solution
    f0_aux <- ECM_star_function(as.matrix(constraints_datasets[[i]][, -1]), constraints_datasets[[i]][, 1])
    ECM_star_aux <- as.numeric(f0_aux[2])  # Extract objective value
    ECM_star <- c(ECM_star, ECM_star_aux)
    
    # Extract OLS beta estimates for constraints and bind row-wise
    Sol_OLS_constraints <- rbind(Sol_OLS_constraints,
                                 as.numeric(f0_aux[3:(ncol(constraints_datasets[[i]][, -1]) + 2)]))
  }
  
  # If no tau vector provided, compute tau_min and tau_max to create tau grid
  if (is.null(vector_tauss)) {
    tau_min <- tau_min_function(constraints_datasets, ECM_star)
    tau_max_lambda_vector <- tau_max(constraints_datasets, vector_lambdas, alpha)
    
    # Create tau sequence from tau_min up to max tau_max + 2 with defined resolution
    vector_tauss <- seq(tau_min, max(tau_max_lambda_vector) + 2, by = grosor_rejilla_taus)
  }
  
  # Initialize counters and lists to store results
  j = 0
  valorobjetivo <- list()  # Store objective values for each parameter combo
  betas <- list()          # Store beta coefficient solutions
  status_vector <- list()  # Store solver statuses
  
  # Loop over alpha values (only one in this case)
  for (alpha in vector_alphas) {
    # Loop over lambda values in grid
    for (lambda in vector_lambdas) {
      j = j + 1
      
      # Temporary storage for each tau iteration
      valorobjetivo_aux <- c()
      betas_aux <- c()
      status_vector_aux <- c()
      
      # Loop over tau values
      for (tau in vector_tauss) {
        
        ## Construct constraint matrix A for the quadratic program
        # Matrix A_aux combines zeros for intercept and identity matrix for coefficients (except intercept)
        A_aux <- cbind(rep(0, ncol(X) - 1), diag(ncol(X) - 1))
        
        # Full constraint matrix includes slack variables (-I and I matrices)
        A <- as.matrix(cbind(A_aux, -diag(nrow(A_aux)), diag(nrow(A_aux))))
        
        # Right-hand side for constraint
        rhs <- 0
        
        # Constraint sense: equality
        sense <- '='
        
        # Constant term in the objective function (scaled squared norm of y)
        objcon <- (1 / nrow(X)) * t(y) %*% y
        
        # Quadratic term matrix Q including ridge penalty weighted by lambda*(1-alpha)
        matrizcompuesta <- (1 / nrow(X)) * as.matrix(
          rbind(
            cbind(t(X) %*% X + lambda * (1 - alpha) * (t(A_aux) %*% A_aux),
                  matrix(rep(0, ncol(X) * 2 * nrow(A)), ncol = 2 * nrow(A))),
            matrix(rep(0, 2 * nrow(A) * (ncol(X) + 2 * nrow(A))),
                   nrow = 2 * nrow(A))
          )
        )
        Q <- matrizcompuesta
        
        # Linear term vector in the objective function:
        # Includes L1 penalty alpha*lambda for slack variables and negative gradient term from least squares
        obj <- c(rep(0, ncol(X)), rep(alpha * lambda, 2 * nrow(A))) -
          (2 / nrow(X)) * t(y) %*% cbind(X, matrix(rep(0, 2 * nrow(A) * nrow(X)), ncol = 2 * nrow(A)))
        
        # Variable type: continuous
        vtype <- 'C'
        
        # Minimization problem
        modelsense <- 'min'
        
        # Model name (optional)
        modelname <- 'CLassoMatrizId'
        
        # Lower bounds: no lower bound for beta, zero lower bound for slack variables
        lb <- c(rep(-Inf, ncol(X)), rep(0, 2 * nrow(A)))
        
        # Initialize quadratic constraints list
        quadcon <- list()
        Cuadratica <- list()
        
        # For each constraint dataset, define quadratic constraint representing ECM* condition
        for (i in 1:length(constraints_datasets)) {
          matrizcompuestarestriccion <- as.matrix(
            rbind(
              cbind(t(constraints_datasets[[i]][, -1]) %*% constraints_datasets[[i]][, -1],
                    matrix(rep(0, ncol(constraints_datasets[[i]][, -1]) * 2 * nrow(A)), ncol = 2 * nrow(A))),
              matrix(rep(0, 2 * nrow(A) * (ncol(constraints_datasets[[i]][, -1]) + 2 * nrow(A))),
                     nrow = 2 * nrow(A))
            )
          )
          
          # Quadratic coefficient matrix for constraint i
          Cuadratica$Qc[[i]] <- (1 / nrow(constraints_datasets[[i]][, -1])) * matrizcompuestarestriccion
          
          # Linear coefficient vector for constraint i
          Cuadratica$q[[i]] <- -(2 / nrow(constraints_datasets[[i]][, -1])) *
            t(constraints_datasets[[i]][, 1]) %*% cbind(constraints_datasets[[i]][, -1],
                                                        matrix(rep(0, 2 * nrow(A) * nrow(constraints_datasets[[i]][, -1])), ncol = 2 * nrow(A)))
          
          # Right hand side for constraint i, includes tau relaxation term
          Cuadratica$rhs[[i]] <- (1 + tau) * ECM_star[i] -
            (1 / nrow(constraints_datasets[[i]][, -1])) * t(constraints_datasets[[i]][, 1]) %*% constraints_datasets[[i]][, 1]
          
          # Inequality sense (less than or equal)
          Cuadratica$sense[[i]] <- '<'
          
          # Append quadratic constraint i to the list in Gurobi format
          quadcon[[i]] <- list(Qc = Cuadratica$Qc[[i]],
                               q = as.numeric(Cuadratica$q[[i]]),
                               rhs = Cuadratica$rhs[[i]],
                               sense = Cuadratica$sense[[i]])
        }
        
        # Build full model list for Gurobi solver including constraints and quadratic constraints
        model <- list(A = A,
                      rhs = rhs,
                      sense = sense,
                      obj = obj,
                      Q = Q,
                      vtype = vtype,
                      modelsense = modelsense,
                      lb = lb,
                      objcon = objcon,
                      quadcon = quadcon)
        
        # Solve the quadratic programming problem with Gurobi
        result <- gurobi(model)
        
        # Store optimization status
        status_vector_aux <- c(status_vector_aux, result$status)
        
        # Extract objective value and beta coefficients (first ncol(X) variables)
        valorobjetivo_aux1 <- result$objval
        betas_aux1 <- result$x[1:ncol(X)]
        
        # Append results for current tau
        valorobjetivo_aux <- c(valorobjetivo_aux, valorobjetivo_aux1)
        betas_aux <- rbind(betas_aux, betas_aux1)
      }
      
      # Store results for current (alpha, lambda) combination
      valorobjetivo[[j]] <- valorobjetivo_aux
      betas[[j]] <- betas_aux
      status_vector[[j]] <- status_vector_aux
    }
  }
  
  # Return a list with beta solutions and tau grid used
  return(list(betas = betas, taus = vector_tauss))
}
