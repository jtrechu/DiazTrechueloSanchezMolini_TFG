tableData <- function(){
  # Set number of selected covariates and total covariates
  K <- 8
  p <- 20
  
  # Define a p x p covariance matrix with correlation decaying exponentially with distance
  Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
  
  # Define sample sizes: 3 groups with 100 samples each, and 5 groups with 500 samples each
  n <- c(rep(400,3),rep(1200,5))
  N <- sum(n)  # Total number of samples
  
  # Define two different beta coefficient vectors for different sets of observations
  betas1 <- c(rep(1,3), rep(1+8^.5,5))      # Used for first 300 samples
  betas2 <- c(rep(1+8^.5,3), rep(1,5))      # Used for remaining 2500 samples
  
  # Load the MASS package for multivariate normal generation
  library(MASS)
  set.seed(2025)  # Set seed for reproducibility
  
  # Generate N samples from a multivariate normal distribution with mean 0 and covariance Sigma
  samples <- mvrnorm(n = N, mu = rep(0,20) , Sigma = Sigma)
  
  # Randomly select K (8) covariates to use in generating response variable
  covariates <- sample(1:20, 8)
  
  # Initialize response vector Y
  Y <- c()
  
  # Generate response values for first 300 samples using betas1 and additive Gaussian noise
  for (i in 1:1200) {
    Y <- c(Y,sum(samples[i,covariates]*betas1)+rnorm(1))
  }
  
  # Generate response values for the remaining 2500 samples using betas2 and additive Gaussian noise
  for (i in 1201:7200) {
    Y <- c(Y,sum(samples[i,covariates]*betas2)+rnorm(1))
  }
  
  # Convert samples to a data frame
  samples <- as.data.frame(samples)
  
  # Add response variable Y to the data frame
  samples$Y <- Y
  
  # Add group labels: Groups 1 to 3 have 100 samples each, groups 4 to 8 have 500 samples each
  samples$Group <- as.factor(c(rep(1:3,each=400),rep(4:8,each=1200)))
  
  # Add an intercept column to the data frame
  samples <- cbind(Intercept=1, samples)
  
  # Return the complete data frame
  return(samples)
}
