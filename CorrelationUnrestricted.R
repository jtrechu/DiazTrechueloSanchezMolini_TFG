source("CorrelationDataSet.R")
ext_lambdas <- 30

grosor_rejilla_lambdas = 0.1
vector_lambdas <- seq(0, ext_lambdas, by = grosor_rejilla_lambdas)
df <- generate_dataset(0.9, 100, 80)
df <- df[, -12]

library(glmnet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(latex2exp)

# Prepare data
x <- as.matrix(df[, -1])
y <- df[, 1]

# Fit Elastic Net model
fit <- glmnet(
  x, y,
  alpha = 1,
  lambda = seq(10, 0, length.out = 100),
  standardize = FALSE
)

coefs <- as.matrix(fit$beta)
lambda_vals <- fit$lambda

# Build data frame
coef_df <- as.data.frame(t(coefs)) %>%
  mutate(lambda = lambda_vals) %>%
  pivot_longer(
    cols = -lambda,
    names_to = "variable",
    values_to = "coefficient"
  )

vars <- unique(coef_df$variable)
vars_ordered <- c(setdiff(vars, "X10"), "X10")

coef_df$variable <- factor(coef_df$variable, levels = vars_ordered)

ggplot(coef_df, aes(x = lambda, y = coefficient, color = variable)) +
  geom_line(size = 1.2) +  # thicker lines
  theme_minimal() +
  labs(
    title = TeX("Regularization Paths for the Lasso ($\\alpha=1$)"),
    x = expression(lambda),
    y = TeX("$\\beta_i$")
  )
