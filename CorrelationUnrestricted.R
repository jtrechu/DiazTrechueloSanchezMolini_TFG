# Load the dataset generation function
source("CorrelationDataSet.R")

# Set maximum lambda value and grid step size
ext_lambdas <- 30
grosor_rejilla_lambdas <- 0.1  # Grid step size for lambda

# Create sequence of lambda values from 0 to 30 with step 0.1
vector_lambdas <- seq(0, ext_lambdas, by = grosor_rejilla_lambdas)

# Generate dataset with:
# - Correlation p = 0.9 between X1 and X2
# - 100 total observations
# - 80 in group 1, 20 in group 2
df <- generate_dataset(0.9, 100, 80)

# Remove the type column (12th column) since we're doing standard EN
df <- df[, -12]

# Load required libraries
library(glmnet)    # For Elastic Net implementation
library(ggplot2)   # For plotting
library(dplyr)     # For data manipulation
library(tidyr)     # For data reshaping
library(latex2exp) # For LaTeX expressions in plots

# Prepare data matrices:
# - x: predictor matrix (all columns except first)
# - y: response variable (first column)
x <- as.matrix(df[, -1])
y <- df[, 1]

# Fit Lasso model (alpha=1 is pure Lasso)
# Note: We're using decreasing lambda sequence for better path visualization
fit <- glmnet(
  x, y,
  alpha = 1,               # 1 = Lasso, 0 = Ridge, (0,1) = Elastic Net
  lambda = seq(10, 0, length.out = 100),  # Decreasing sequence for path
  standardize = FALSE       # Don't standardize (since we want raw coefficients)
)

# Extract coefficients (excluding intercept)
coefs <- as.matrix(fit$beta)

# Get corresponding lambda values
lambda_vals <- fit$lambda

# Transform coefficients into tidy data frame for ggplot:
# 1. Transpose coefficients matrix
# 2. Convert to data frame and add lambda column
# 3. Reshape to long format (one row per coefficient per lambda)
coef_df <- as.data.frame(t(coefs)) %>%
  mutate(lambda = lambda_vals) %>%
  pivot_longer(
    cols = -lambda,
    names_to = "variable",
    values_to = "coefficient"
  )

# Get unique variable names and reorder to put X10 last
vars <- unique(coef_df$variable)
vars_ordered <- c(setdiff(vars, "X10"), "X10")  # All vars except X10, then X10

# Convert variable to factor with specified ordering
coef_df$variable <- factor(coef_df$variable, levels = vars_ordered)

# Create coefficient path plot
ggplot(coef_df, aes(x = lambda, y = coefficient, color = variable)) +
  geom_line(size = 1.2) +  # Thicker lines for better visibility
  theme_minimal() +        # Clean theme
  labs(
    title = TeX("Regularization Paths for the Lasso ($\\alpha=1$)"),
    x = expression(lambda),  # Proper lambda symbol
    y = TeX("$\\beta_i$")   # Proper beta symbol with subscript
  ) 
