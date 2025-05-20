# Load required packages
library(glmnet)    # For fitting Lasso and Elastic Net models
library(ggplot2)   # For creating visualizations
library(reshape2)  # For data manipulation (though not directly used here)

# Set seed for reproducibility of random data generation
set.seed(2025)

# ==============================================
# SIMULATION PARAMETERS
# ==============================================
# We create a high-dimensional setting where p > n
# This demonstrates the saturation property of Lasso
n <- 10           # Number of observations (small on purpose)
p <- 20            # Number of predictors (larger than n)
significant_vars <- 15# Number of truly important variables

# ==============================================
# DATA GENERATION
# ==============================================
# Generate predictor matrix with standard normal distribution
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", 1:p)  # Name variables X1 to X20

# Create true coefficients:
# - First 15 variables have β=1 (important)
# - Last 5 variables have β=0 (noise)
beta <- c(rep(1, significant_vars), rep(0, p - significant_vars))

# Generate response variable with noise
epsilon <- rnorm(n)  # Random noise
Y <- X %*% beta + epsilon  # True relationship plus noise

# ==============================================
# MODEL FITTING
# ==============================================
# Fit Lasso model (α=1 for pure Lasso)
# Note: glmnet automatically generates a sensible lambda sequence
lasso_fit <- glmnet(X, Y, alpha = 1)

# Fit Elastic Net model with α=0.2 (20% Lasso penalty, 80% Ridge)
en_fit <- glmnet(X, Y, alpha = 0.2, lambda = lasso_fit$lambda)

# ==============================================
# ANALYSIS OF SELECTED VARIABLES
# ==============================================
# Function to count non-zero coefficients (excluding intercept)
# Threshold accounts for numerical precision in coefficient estimates
count_nonzero <- function(fit, threshold = 1e-8) {
  apply(coef(fit), 2, function(x) sum(abs(x) > threshold)) - 1  # subtract intercept
}

# Count selected variables across all lambda values
lasso_counts <- count_nonzero(lasso_fit)  # For Lasso
en_counts <- count_nonzero(en_fit)        # For Elastic Net

# ==============================================
# PREPARE DATA FOR VISUALIZATION
# ==============================================
plot_data <- data.frame(
  lambda = c(lasso_fit$lambda, en_fit$lambda),
  count = c(lasso_counts, en_counts),
  model = rep(c("Lasso", "Elastic Net (α=0.2)"), 
              c(length(lasso_fit$lambda), length(en_fit$lambda)))
)

# ==============================================
# CREATE THE PLOT
# ==============================================
ggplot(plot_data, aes(x = lambda, y = count, color = model)) +
  geom_line(size = 0.8) +  # Slightly thinner lines for clarity
  # Add reference line at y=n to show saturation point
  geom_hline(yintercept = n, linetype = "dashed", color = "black") +
  # Use log scale for lambda (common for regularization path plots)
  scale_x_log10() +  
  labs(
    x = "Lambda (log scale)",
    y = "Number of selected variables",
    title = "Variable Selection with Lasso and Elastic Net",
    subtitle = paste("Demonstrating Lasso's saturation at n =", n, "variables when p > n")
  ) +
  theme_minimal() +  # Clean, modern theme
  # Custom color scheme for clear differentiation
  scale_color_manual(
    values = c("Lasso" = "blue", "Elastic Net (α=0.2)" = "red"),
    name = "Model Type"
  ) +
  # Additional plot customization
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40")
  )