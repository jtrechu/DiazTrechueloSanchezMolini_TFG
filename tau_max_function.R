tau_max <- function(data_list, lambda_seq,alpha) {
  library(glmnet)
  L <- length(data_list)
  tau_values <- numeric(length(lambda_seq))
  
  for (i in seq_along(lambda_seq)) {
    lambda <- lambda_seq[i]
    ratios <- numeric(L)
    
    for (l in 1:L) {
      data <- data_list[[l]]
      
      
      x <- as.matrix(data[, -c(1,2)])
      y <- data[, 1]
      
      # OLS
      ols_fit <- lm(y ~ x)
      y_pred_ols <- predict(ols_fit)
      mse_ols <- mean((y - y_pred_ols)^2)
      
      # Lasso
      en_fit <- glmnet(x, y, alpha = alpha, lambda = lambda/2, standardize = FALSE) 
      #glmnet minimizes 1/2*RSS/nobs
      pred_en <- predict(en_fit, s = lambda, newx = x)
      mse_en <- mean((y - pred_en)^2)
      
      # Ratio
      ratios[l] <- (mse_en / mse_ols) - 1
    }
    
    tau_values[i] <- max(ratios)
  }
  
  return(tau_values)
}
