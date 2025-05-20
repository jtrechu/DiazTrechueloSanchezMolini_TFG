library(ggplot2)
library(latex2exp)
source('ProstateCSCEN.R')
alpha <- 0
l2 <- CSCEN(lambda = 150, alpha = alpha, grosor_lambda=0.5)
x <- l2$taus 
y <- seq(0,150,by=0.5)
grid <- expand.grid(x = x, y = y)
z<- c()
for (i in 1:length(y)){
  z <- c(z,l2$betas[[i]][,2])
}

tau_min_val <- min(l2$taus)
tau_max_val <- max(l2$taus)-2

ggplot(grid, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "navy",    # valores negativos
    mid = "white",   # cerca de cero
    high = "firebrick",    # valores positivos
    midpoint = 0,    # el valor central blanco
  )+
  scale_x_continuous(
    breaks = c(tau_min_val, tau_max_val),
    labels = c(TeX("$\\tau_{min}$"), TeX("$\\tau_{max}$"))
  ) +
  labs(
    x = TeX("$\\tau$"),
    y = TeX("$\\lambda$"),
    fill = TeX("$\\beta_1^{CSCEN}$")
  ) +
  theme_minimal() +
  ggtitle(TeX("Evolution of $\\beta_1$ for Cost Sensitive Constrained Ridge Regression")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
