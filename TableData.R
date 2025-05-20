##Datos sinteticos 2
tableData <- function(){
  K <- 8
  p <- 20
  Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
  n <- c(rep(100,3),rep(500,5))
  N <- sum(n)
  betas1 <- c(rep(1,3), rep(1+10^.5,5))
  betas2 <- c(rep(1+10^.5,3), rep(1,5))
  
  library(MASS)
  set.seed(2025)
  samples <- mvrnorm(n = N, mu = rep(0,20) , Sigma = Sigma)
  samples
  covariates <- sample(1:20, 8)
  Y <- c()
  for (i in 1:300) {
    Y <- c(Y,sum(samples[i,covariates]*betas1)+rnorm(1))
  }
  for (i in 301:2800) {
    Y <- c(Y,sum(samples[i,covariates]*betas2)+rnorm(1))
  }
  samples <- as.data.frame(samples)
  samples$Y <- Y
  samples$Group <- as.factor(c(rep(1:3,each=100),rep(4:8,each=500)))
  samples <- cbind(Intercept=1, samples)
  return(samples)
}
