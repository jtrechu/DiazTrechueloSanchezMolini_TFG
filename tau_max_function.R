# Function to calculate maximum tau values across different lambda values
# Arguments:
#   data_list: List of datasets
#   lambda_seq: Sequence of lambda values to test
#   alpha: Elastic net mixing parameter (1 for lasso, 0 for ridge)
tau_max <- function(data_list, lambda_seq, alpha) {
  # Load required library for elastic net implementation
  library(glmnet)
  
  # Get number of datasets in the list
  L <- length(data_list)
  
  # Initialize vector to store tau values for each lambda
  tau_values <- numeric(length(lambda_seq))
  
  # Loop through each lambda value in the sequence
  for (i in seq_along(lambda_seq)) {
    lambda <- lambda_seq[i]
    
    # Initialize vector to store ratios for each dataset
    ratios <- numeric(L)
    
    # Loop through each dataset in the list
    for (l in 1:L) {
      data <- data_list[[l]]
      
      # Prepare data: x is all columns except first two, y is first column
      x <- as.matrix(data[, -c(1,2)])
      y <- data[, 1]
      
      # Fit Ordinary Least Squares (OLS) model
      ols_fit <- lm(y ~ x)
      # Predict using OLS model
      y_pred_ols <- predict(ols_fit)
      # Calculate Mean Squared Error for OLS
      mse_ols <- mean((y - y_pred_ols)^2)
      
      # Fit elastic net model (lasso when alpha=1)
      # Note: glmnet uses lambda/2 because it minimizes (1/2)*RSS/nobs + lambda*penalty
      en_fit <- glmnet(x, y, alpha = alpha, lambda = lambda/2, standardize = FALSE) 
      # Predict using elastic net model
      pred_en <- predict(en_fit, s = lambda, newx = x)
      # Calculate Mean Squared Error for elastic net
      mse_en <- mean((y - pred_en)^2)
      
      # Calculate performance ratio: (MSE_en/MSE_ols) - 1
      # This measures relative performance degradation compared to OLS
      ratios[l] <- (mse_en / mse_ols) - 1
    }
    
    # Store the maximum ratio across all datasets for this lambda value
    tau_values[i] <- max(ratios)
  }
  
  # Return vector of maximum tau values for each lambda
  return(tau_values)
}