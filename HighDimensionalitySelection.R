# Load required packages
library(glmnet)
library(ggplot2)
library(reshape2)

# Set seed for reproducibility
set.seed(2025)

# Parameters
n <- 10  # Number of observations
p <- 20   # Number of predictors
significant_vars <- 15  # Number of significant variables

# Generate data
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", 1:p)

# Create response variable (first 15 variables are significant)
beta <- c(rep(1, significant_vars), rep(0, p - significant_vars))
epsilon <- rnorm(n)
Y <- X %*% beta + epsilon

# Fit Lasso model (alpha = 1)
lasso_fit <- glmnet(X, Y, alpha = 1)

# Fit Elastic Net model with alpha = 0.2
en_fit <- glmnet(X, Y, alpha = 0.2)

# Function to count non-zero coefficients for a given lambda
count_nonzero <- function(fit, threshold = 1e-8) {
  apply(coef(fit), 2, function(x) sum(abs(x) > threshold)) - 1  # subtract intercept
}

# Get number of selected variables for each lambda
lasso_counts <- count_nonzero(lasso_fit)
en_counts <- count_nonzero(en_fit)

# Create data frame for plotting
plot_data <- data.frame(
  lambda = c(lasso_fit$lambda, en_fit$lambda),
  count = c(lasso_counts, en_counts),
  model = rep(c("Lasso", "Elastic Net (α=0.2)"), 
              c(length(lasso_fit$lambda), length(en_fit$lambda)))
)

# Create plot
ggplot(plot_data, aes(x = lambda, y = count, color = model)) +
  geom_line(size = 0.8) +  # thinner lines
  geom_hline(yintercept = n, linetype = "dashed", color = "black") +
  scale_x_log10() +  # Log scale for lambda
  labs(
    x = "Lambda (log scale)",
    y = "Number of selected variables",
    title = "Variable Selection with Lasso and Elastic Net"
  ) +
  theme_minimal() +  # white background
  scale_color_manual(values = c("Lasso" = "blue", "Elastic Net (α=0.2)" = "red"))
