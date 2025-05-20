# Load required libraries
library(ggplot2)      # For creating plots
library(latex2exp)    # For LaTeX mathematical notation in plots

# Load custom CSCEN function
source('ProstateCSCEN.R')

# Set alpha parameter (0 = pure ridge regression, no lasso component)
alpha <- 0

# Run CSCEN function with:
# - lambda up to 150
# - specified alpha
# - lambda step size of 0.5
l2 <- CSCEN(lambda = 150, alpha = alpha, grosor_lambda = 0.5)

# Extract tau values from CSCEN results
x <- l2$taus 

# Create sequence of lambda values from 0 to 150 in steps of 0.5
y <- seq(0, 150, by = 0.5)

# Create a grid of all combinations of tau and lambda values
grid <- expand.grid(x = x, y = y)

# Extract the second coefficient (lcavol) for all lambda values
z <- c()
for (i in 1:length(y)) {
  z <- c(z, l2$betas[[i]][, 2])  # [,2] selects the second coefficient (lcavol)
}

# Determine tau range for x-axis labels
tau_min_val <- min(l2$taus)       # Minimum tau value
tau_max_val <- max(l2$taus) - 2   # Maximum tau value minus 2 (for better visualization)

# Create the heatmap plot
ggplot(grid, aes(x = x, y = y, fill = z)) +
  # Create tile plot (heatmap)
  geom_tile() +
  
  # Define color gradient for the heatmap:
  # - navy for negative values
  # - white near zero
  # - firebrick for positive values
  scale_fill_gradient2(
    low = "navy",     # Color for lowest values
    mid = "white",    # Color for midpoint values
    high = "firebrick",  # Color for highest values
    midpoint = 0      # Value that should be colored with 'mid' color
  ) +
  
  # Customize x-axis:
  # - Only show min and max tau values
  # - Label them with LaTeX symbols
  scale_x_continuous(
    breaks = c(tau_min_val, tau_max_val),
    labels = c(TeX("$\\tau_{min}$"), TeX("$\\tau_{max}$"))
  ) +
  
  # Set axis labels and legend title using LaTeX notation
  labs(
    x = TeX("$\\tau$"),                  # x-axis label
    y = TeX("$\\lambda$"),               # y-axis label
    fill = TeX("$\\beta_1^{CSCEN}$")     # Legend title (lcavol coefficient)
  ) +
  
  # Use minimal theme for clean appearance
  theme_minimal() +
  
  # Add plot title with LaTeX notation
  ggtitle(TeX("Evolution of $\\beta_1$ for Cost Sensitive Constrained Ridge Regression")) +
  
  # Adjust x-axis text appearance
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Keep labels horizontal and centered