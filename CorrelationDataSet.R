# Load the MASS package for multivariate normal distribution generation
library(MASS)  # Provides mvrnorm() function

# Function to generate a synthetic dataset with correlated predictors
# Parameters:
#   p: Correlation coefficient between X1 and X2
#   n: Total number of observations
#   n1: Number of observations in group 1 (remainder will be in group 2)
generate_dataset <- function(p, n, n1) {
  # Set random seed for reproducibility
  set.seed(2025)
  
  # Calculate number of observations in group 2
  n2 <- n - n1
  n_total <- n  # Total observations (same as input n, but for clarity)
  
  # Create covariance matrix for X1 and X2:
  # - Variance of each variable = 1 (diagonal)
  # - Covariance between X1 and X2 = p (off-diagonal)
  Sigma <- matrix(c(
    1, p,  # First row: X1 variance, X1-X2 covariance
    p, 1   # Second row: X1-X2 covariance, X2 variance
  ), nrow = 2, byrow = TRUE)
  
  # Generate X1 and X2 from multivariate normal distribution with:
  # - n_total observations
  # - Mean vector of zeros
  # - Covariance matrix Sigma
  X12 <- mvrnorm(n = n_total, mu = rep(0, 2), Sigma = Sigma)
  
  # Generate 8 additional independent predictors (X3 to X10) from standard normal
  X3_10 <- matrix(rnorm(n_total * 8), ncol = 8)
  
  # Combine all predictors into one matrix
  X <- cbind(X12, X3_10)
  # Name the columns X1 through X10
  colnames(X) <- paste0("X", 1:10)
  
  # Generate random noise/error term from standard normal
  eps <- rnorm(n, mean = 0, sd = 1)
  
  # Create response variable Y as a linear combination of predictors:
  # Y = X1 + X2 - 3*X3 + error
  # Note: This creates a true underlying relationship where:
  # - X1 and X2 have positive effects
  # - X3 has a strong negative effect
  # - X4-X10 are irrelevant (not included in the true model)
  Y <- X[,1] + X[,2] - 3*X[,3] + eps
  
  # Create grouping variable:
  # - First n1 observations are type 1
  # - Remaining n2 observations are type 2
  type <- c(rep(1, n1), rep(2, n2))
  
  # Combine into a data frame with:
  # - Y (response)
  # - X1-X10 (predictors)
  # - type (grouping factor)
  df <- data.frame(Y, X, type = as.factor(type))
  
  return(df)
}