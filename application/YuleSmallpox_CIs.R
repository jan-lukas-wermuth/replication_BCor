# ===========================================================
# Load packages and source files
# ===========================================================
library(epitools)
library(dplyr)
library(DescTools)
library(BCor)

# Calculate confidence intervals for Introductory example (Yule 1912)

# Create table
Smallpox_table <- as.table(matrix(c(197,139,2,19), nrow=2))
rownames(Smallpox_table) <- c(1,0) 
colnames(Smallpox_table) <- c(1,0)

Smallpox_data <- expand.table(Smallpox_table)
X <- as.numeric(Smallpox_data[,1])-1
Y <- as.numeric(Smallpox_data[,2])-1

BCor::YuleQ(X, Y, alpha = 0.9, Fisher = FALSE)
BCor::YuleQ(X, Y, alpha = 0.9, Fisher = TRUE)
BCor::Cole(X, Y, alpha = 0.9, Fisher = FALSE)
BCor::Cole(X, Y, alpha = 0.9, Fisher = TRUE)
BCor::Phi(X, Y, alpha = 0.9, Fisher = FALSE)
BCor::Phi(X, Y, alpha = 0.9, Fisher = TRUE)

# Results: 
# Q_est = 0.862 (0.702, 1.02) vs. Z: (0.592, 0.958)
# C_est = 0.829 (0.606, 0.999) vs. Z: (0.438, 0.956)
# Phi_est = 0.233 (0.164, 0.301) vs. Z: (0.163, 0.300)

DescTools::CramerV(Smallpox_table) # 0.23
DescTools::ContCoef(Smallpox_table) # 0.23
BCor::Tetrachoric(Smallpox_table) # 0.61
DescTools::YuleY(Smallpox_table) # 0.57
