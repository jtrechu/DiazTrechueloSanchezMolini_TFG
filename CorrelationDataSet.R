library(MASS)  # Para mvrnorm

generate_dataset <- function(p,n,n1) {
  set.seed(2025)
  
  # Número de observaciones
  n2 <- n-n1
  n_total <- n
  
  # Matriz de covarianzas para X1, X2, X3
  Sigma <- matrix(c(
    1, p,
    p, 1
  ), nrow = 2, byrow = TRUE)
  
  # Generar X1, X2, X3 ~ N(0, Sigma)
  X12 <- mvrnorm(n = n_total, mu = rep(0, 2), Sigma = Sigma)
  
  # Generar X4 a X20 ~ N(0, 1) independientes
  X3_10 <- matrix(rnorm(n_total * 8), ncol = 8)
  # Combinar todas las X
  X <- cbind(X12, X3_10)
  colnames(X) <- paste0("X", 1:10)
  
  # Generar el error
  eps <- c(rnorm(n, mean = 0, sd = 1))
  
  # Calcular Y = X1 + X2 + X3 - 3*X4 + error
  Y <- X[,1] + X[,2] -3*X[,3] + eps
  
  
  # Crear variable de tipo
  type <- c(rep(1, n1), rep(2, n2))
  
  # Crear dataframe
  df <- data.frame(Y, X, type = as.factor(type))
  return(df)
}



