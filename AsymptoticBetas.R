
#We'll take TAU = 1.5  and check that beta=0 is not feasible
library(genridge)
data(prostate)
Prostate<-cbind(Intercept = 1, prostate[,-10])
X <- Prostate[,-10]
X <- as.matrix(X)
y <-  Prostate[,10]
X_menoresde65 <- Prostate[which(Prostate[,"age"]<65),-10]
X_mayoresde65 <- Prostate[which(Prostate[,"age"]>=65),-10]
y_menoresde65 <- Prostate[which(Prostate[,"age"]<65),10]
y_mayoresde65 <- Prostate[which(Prostate[,"age"]>=65),10]

ErrorIfBeta0 <- (1/length(y_menoresde65))*sum((y_menoresde65)^2)

model <- lm(y_menoresde65 ~ .,data = X_menoresde65[,-1])  
rhs <- 2.5*(mean(residuals(model)^2))
ErrorIfBeta0<=rhs ##NOT FEASIBLE


source("ProstateCSCEN.R")
l <- CSCEN(lambda = 100,alpha = 0,grosor_lambda=0.5,vector_tauss = c(1.5))
l$taus

library(ggplot2)
library(latex2exp)
library(tidyr)
library(dplyr)

# Convert list to data frame
lambda <- seq(0,100,by=0.5)
betas <- l$betas
beta_matrix <- do.call(rbind, betas)
beta_matrix <- cbind(lambda,beta_matrix[,-1])
colnames(beta_matrix) <- c("lambda", "lcavol","lweight","age","lbph","svi","lcp","gleason","pgg45")
df <- as.data.frame(beta_matrix)

df_long <- pivot_longer(df, cols = c("lcavol","lweight","age","lbph","svi","lcp","gleason","pgg45"), 
                        names_to = "coefficient", values_to = "value")

# Plot
ggplot(df_long, aes(x = lambda, y = value, color = coefficient)) +
  geom_line(linewidth = 1) +
  labs(x = TeX("$\\lambda$"), y = TeX("$\\beta_j$")) +
  theme_minimal() +
  ggtitle(TeX("Evolution of $\\underline{\\beta}^{CSCRidge}(\\lambda)$ with $\\tau=1.5$"))+
  theme(legend.title = element_blank())


