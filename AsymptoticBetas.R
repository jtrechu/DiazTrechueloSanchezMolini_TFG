# Load required libraries and data
library(genridge)
data(prostate)

# Prepare the Prostate dataset - add intercept column and separate predictors from response
Prostate <- cbind(Intercept = 1, prostate[,-10])  # Add column of 1s for intercept
X <- Prostate[,-10]  # Predictor matrix
X <- as.matrix(X)    # Convert to matrix format
y <- Prostate[,10]   # Response variable

# Split data into two age groups: under 65 and 65+
X_menoresde65 <- Prostate[which(Prostate[,"age"]<65),-10]  # Predictors for age < 65
X_mayoresde65 <- Prostate[which(Prostate[,"age"]>=65),-10] # Predictors for age >= 65
y_menoresde65 <- Prostate[which(Prostate[,"age"]<65),10]   # Response for age < 65
y_mayoresde65 <- Prostate[which(Prostate[,"age"]>=65),10]  # Response for age >= 65

# Calculate the error when using beta=0 (null model) for younger age group
ErrorIfBeta0 <- (1/length(y_menoresde65))*sum((y_menoresde65)^2)

# Fit a linear model to the younger age group data (excluding intercept column)
model <- lm(y_menoresde65 ~ ., data = X_menoresde65[,-1])  

# Calculate the right-hand side of the feasibility condition
rhs <- 2.5*(mean(residuals(model)^2))

# Check if beta=0 is feasible (should NOT be feasible for our analysis)
ErrorIfBeta0 <= rhs  # This should return FALSE

# Load our custom CSCEN function
source("ProstateCSCEN.R")

# Call CSCEN function with specific parameters:
# lambda = 100 (maximum lambda value)
# alpha to choose
# grosor_lambda = 0.5 (step size for lambda sequence)
# vector_tauss = c(1.5) (tau value to use)
l <- CSCEN(lambda = 100, alpha = 0, grosor_lambda = 0.5, vector_tauss = c(1.5))
l$taus  # Display the tau values used (should just show 1.5)

# Load libraries for plotting and data manipulation
library(ggplot2)
library(latex2exp)  # For LaTeX-style mathematical notation in plots
library(tidyr)      # For data reshaping
library(dplyr)      # For data manipulation

# Prepare data for plotting:
# Create sequence of lambda values from 0 to 100 in steps of 0.5
lambda <- seq(0, 100, by = 0.5)

# Extract the beta coefficients from the CSCEN results
betas <- l$betas

# Combine all beta matrices into one matrix and add lambda column
beta_matrix <- do.call(rbind, betas)  # Stack all beta matrices vertically
beta_matrix <- cbind(lambda, beta_matrix[,-1])  # Add lambda column, remove intercept

# Name the columns appropriately
colnames(beta_matrix) <- c("lambda", "lcavol", "lweight", "age", "lbph", 
                           "svi", "lcp", "gleason", "pgg45")

# Convert to data frame for plotting
df <- as.data.frame(beta_matrix)

# Reshape data from wide to long format for ggplot2
# Each coefficient gets its own row with lambda, coefficient name, and value
df_long <- pivot_longer(df, 
                        cols = c("lcavol", "lweight", "age", "lbph", 
                                 "svi", "lcp", "gleason", "pgg45"), 
                        names_to = "coefficient", 
                        values_to = "value")

# Create the plot showing evolution of coefficients with lambda
ggplot(df_long, aes(x = lambda, y = value, color = coefficient)) +
  geom_line(linewidth = 1) +  # Draw lines for each coefficient
  labs(x = TeX("$\\lambda$"), y = TeX("$\\beta_j$")) +  # Axis labels with LaTeX
  theme_minimal() +  # Clean, minimal theme
  ggtitle(TeX("Evolution of $\\underline{\\beta}^{CSCRidge}(\\lambda)$ with $\\tau=1.5$")) +
  theme(legend.title = element_blank())  # Remove legend title
