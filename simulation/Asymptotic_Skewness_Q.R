# ===========================================================
# Load packages and source files
# ===========================================================
library(MASS)
library(DescTools)
library(EnvStats)
library(pracma)
library(dplyr)
library(ggplot2)
library(BCor)

MC <- 100000
n <- 100

Q <- 0.95
QZ <- atanh(Q)
p <- 0.5
q <- 0.5
r <- YuleQ.inv(Q, c(p, 1 - p, q, 1 - q))[1,1]

Q_vec <- rep(NA, MC)
QZ_vec <- rep(NA, MC)

tbl_coding <- tibble(X=c(1,0,1,0), Y=c(1,1,0,0))

for (i in 1:MC) {
  # Simulate
  set.seed(i)
  prob <- c(r, q-r, p-r, 1-p-q+r) %>% pmax(0) %>% pmin(1)
  
  sample_index <- sample(1:4, size=n, prob=prob,  replace = TRUE)
  sample_tbl <- tbl_coding[sample_index,]

  Q_vec[i] <- BCor::YuleQ(sample_tbl$X, sample_tbl$Y, alpha = FALSE)$Q
  QZ_vec[i] <- atanh(Q_vec[i])
}

# Asymptotic Distribution of Q_hat
Omega <- matrix(c(p*(1-p), r - p * q, r * (1-p), r - p * q, q*(1-q), r*(1-q), r*(1-p), r*(1-q), r*(1-r)), ncol = 3)
G_mat <- G(p, q, r)
H_mat <- H(p, q, r)
Q_Var <- as.numeric(t(G_mat) %*% Omega %*% G_mat) / n
QZ_Var <- as.numeric(t(H_mat) %*% Omega %*% H_mat) / n

options(ggplot2.discrete.colour= c("blue", "red", "black"))
ggplot(data.frame(x = 1:100), aes(x = x)) +
  stat_function(aes(colour = "Standard"), fun = function(x) dnorm(x, mean = Q, sd = sqrt(Q_Var)), linewidth = 1) + 
  stat_function(aes(colour = "Fisher"), fun = function(x) {dnorm(atanh(x), mean = QZ, sd = sqrt(QZ_Var)) * abs(1/(1-x^2))}, linewidth = 1) + 
  geom_vline(aes(xintercept = Q, colour = "True Q"), linewidth = 1, show.legend = FALSE) +
  scale_colour_manual(values = c("red", "blue", "black")) +
  geom_histogram(data = data.frame(Q_vec), aes(x = Q_vec, y = ..density..), colour = "black", alpha = 0.2, bins = 18, boundary = 1) +
  geom_vline(aes(xintercept = Q + sqrt(Q_Var)*qnorm(0.05)), color = "red", alpha = 0.5, linewidth = 1, linetype = 2) +
  geom_vline(aes(xintercept = Q + sqrt(Q_Var)*qnorm(0.95)), color = "red", alpha = 0.5, linewidth = 1, linetype = 2) +
  geom_vline(aes(xintercept = tanh(QZ + sqrt(QZ_Var)*qnorm(0.05))), color = "blue", alpha = 0.5, linewidth = 1, linetype = 2) +
  geom_vline(aes(xintercept = tanh(QZ + sqrt(QZ_Var)*qnorm(0.95))), color = "blue", alpha = 0.5, linewidth = 1, linetype = 2) +
  xlim(0.8, 1.05) + xlab("Q") + ylab("Density") +
  scale_colour_discrete(name = "") +
  theme_classic(base_size = 15)

ggsave(filename = "~/Dropbox/DimitriadisPohleWermuth/Binary Correlation/replication_BCor/simulation/output/Q_Simulation.pdf", width = 200, height = 150, device = "pdf", units = "mm")





# Helper functions for the Jacobians used for the Delta-Method
# Yule's Q 
G <- function(p, q, r) {
  denom <- (p * (q - 2 * r) + r * (1 - 2 * q + 2 * r)) ^ 2
  G1 <- 2 * (q - 1) * r * (q - r) / denom
  G2 <- 2 * (p - 1) * r * (p - r) / denom
  G3 <- -(2 * (q * p ^ 2 + p * q * (-1 + q - 2 * r) + r ^ 2)) / denom
  return(c(G1, G2, G3))
}
# Yule's Q with Fisher transformation
H <- function(p, q, r) {
  denom1 <- -2 * (p - r) * (q - r)
  denom2 <- -2 * (1 - p - q + r)
  H1 <- (q - r) / denom1 + 1 / denom2
  H2 <- (p - r) / denom1 + 1 / denom2
  H3 <- (2 * r - p - q) / denom1 - (1 - p - q + 2 * r) / (r * denom2)
  return(c(H1, H2, H3))
}






