# ===========================================================
# Load packages and source files
# ===========================================================
library(tidyverse)
library(tibble)
library(doParallel)
library(foreach)
library(fdrtool)
library(sandwich)
library(MASS)
library(BCor)


### This files runs on 10 kernels with 10k replications for approximately 1.5h


### Illustrate the asymptotic distributions for these cases:
setting_list <- list(c(0.4, 0.8, -0.5),
                     c(0.4, 0.8, 0),
                     c(0.4, 0.8, 0.5),
                     c(0.2, 0.2, -0.5),
                     c(0.2, 0.2, 0),
                     c(0.2, 0.2, 0.5),
                     c(0.3, 0.7, -0.5),
                     c(0.3, 0.7, 0),
                     c(0.3, 0.7, 0.5))

setting_type <- c("normal_neg", "zero", "normal_pos",
                  "normal_neg", "zero", "equal_pos",
                  "equal_neg", "zero", "normal_pos")

nn <- 2000
M_par <- 100
M_internal <- 100

# Loop in parallel
core.max <- 10
start.time <- Sys.time()
cl <- makeCluster(min(parallel::detectCores()-1, M_par, core.max) )
registerDoParallel(cl)
res_df_MC <- foreach(
  i_MC = 1:M_par,
  .combine=rbind,
  .packages=c("tidyverse", "tibble", "foreach", "BCor"),
  .errorhandling="pass"
)%dopar%{
  set.seed(i_MC) # set seed for reproducibility

  BGov_tbl <- tibble()
  tbl_coding <- tibble(X=c(1,0,1,0), Y=c(1,1,0,0))
  
  for (i_internal in 1:M_internal){
    for (i_set in 1:length(setting_list)){
      p <- setting_list[[i_set]][1]
      q <- setting_list[[i_set]][2]
      C <- setting_list[[i_set]][3]
      type <- setting_type[i_set]
      
      m_plus <- min(p,q) - p*q
      m_minus <- p*q - max(0, p+q-1)
      m_joint <- 1*(C >= 0)*m_plus + 1*(C < 0)*m_minus 
      r <- C*m_joint + p*q
      
      # Simulate
      prob <- c(r, q-r, p-r, 1-p-q+r) %>% pmax(0) %>% pmin(1)
      
      sample_index <- sample(1:4, size=nn, prob=prob,  replace = TRUE)
      sample_tbl <- tbl_coding[sample_index,]
      
      C_est <- Cole(sample_tbl$X, sample_tbl$Y, alpha=alpha, Fisher=FALSE, m_rep = 10000)
      if (!("CI_lower" %in% names(C_est))){ C_est <- C_est %>% mutate(CI_lower=NA, CI_upper=NA) }
      
      # Mutate tibbles to save them with all information
      C_est_mod <- as_tibble(C_est) %>% 
        rename(Estimate=C) %>% 
        mutate(i_MC=i_MC, i_internal=i_internal, r=r, p=p, q=q, n=nn, alpha=alpha, Fisher_Trans=FALSE,
               measure="C", 
               Population=1*(r-p*q >= 0) * (r-p*q)/(min(p,q) - p*q) + 1*(r-p*q < 0) * (r-p*q)/(p*q - max(0, p+q-1)))
      
      # Merge
      BGov_tbl <- bind_rows(BGov_tbl,C_est_mod) 
    }
  }
  BGov_tbl
}  
stopCluster(cl)
end.time <- Sys.time()
(run.time <- end.time-start.time)

head(res_df_MC)
saveRDS(res_df_MC, file = "simulation/data/sim_AsyDists_20240223.rds")







### Analytically compute the asymptotic distributions with true p,q,r values

# Load some functions we need:
Jh_pos <- function(p, q, r) {
  Jh <- cbind(c(-q,-p, 1),
              c(1 * (p < q) - q, 1 * (p > q) - p, 0))
  return(Jh)
}

Jh_neg <- function(p, q, r) {
  Jh <- cbind(c(-q,-p, 1),
              c(q - 1 * (p + q > 1), p - 1 * (p + q > 1), 0))
  return(Jh)
}


df_density <- tibble()
for (i_set in 1:length(setting_list)){
  p <- setting_list[[i_set]][1]
  q <- setting_list[[i_set]][2]
  C <- setting_list[[i_set]][3]
  type <- setting_type[i_set]
  
  m_plus <- min(p,q) - p*q
  m_minus <- p*q - max(0, p+q-1)
  m_joint <- 1*(C >= 0)*m_plus + 1*(C < 0)*m_minus 
  r <- C*m_joint + p*q
  
  # Asymptotic distribution
  Omega <- rbind(c(p*(1-p), r-p*q,   r*(1-p)),
                 c(r-p*q,   q*(1-q), r*(1-q)),
                 c(r*(1-p), r*(1-q), r*(1-r)))
  
  ### Inference for Cole Coefficient C
  # The case {C=0}
  Delta <- c(-q, -p, 1)
  C_Var_zero_pos <-  1/m_plus^2 * as.numeric(t(Delta) %*% Omega %*% Delta)/nn
  C_Var_zero_neg <- 1/m_minus^2 * as.numeric(t(Delta) %*% Omega %*% Delta)/nn
  
  # The cases {C>0 and p!=q} and {C<0 and p!=1-q}
  Gamma_pos <- c(1/m_plus, -(r-p*q)/m_plus^2)
  Gamma_neg <- c(1/m_minus, -(r-p*q)/m_minus^2)
  Jh_pos_est <- Jh_pos(p, q, r)
  Jh_neg_est <- Jh_neg(p, q, r)
  C_Var_pos <- as.numeric(t(Gamma_pos) %*% t(Jh_pos_est) %*% Omega %*% Jh_pos_est %*% Gamma_pos)/nn
  C_Var_neg <- as.numeric(t(Gamma_neg) %*% t(Jh_neg_est) %*% Omega %*% Jh_neg_est %*% Gamma_neg)/nn
  
  
  # The cases {C>0 and p=q} and {C<0 and p=1-q}
  m_rep <- 10^6
  Jf <- rbind(c(1,0,0),
              c(0,1,0),
              c(q,p,0),
              c(-q,-p,1))
  V_Var <- Jf %*% Omega %*% t(Jf)
  
  # Draw from the (estimated) asymptotic distribution
  mV <- MASS::mvrnorm(m_rep, mu=rep(0,4), Sigma=V_Var) # Generate a degenerate normal distribution
  
  sigma_mpos_vec <- rbind(mV[,4],
                          1*(mV[,1] <= mV[,2]) * (mV[,1]-mV[,3]) + 1*(mV[,1] > mV[,2]) * (mV[,2]-mV[,3])) 
  Cdiff_pos_samples <- as.numeric(Gamma_pos %*% sigma_mpos_vec)   # Samples from sqrt(n) (hat C_n - C)

  sigma_mneg_vec <- rbind(mV[,4],
                          1*(mV[,1] <= -mV[,2]) * mV[,3] + 1*(mV[,1] > -mV[,2]) * (mV[,3] - mV[,1] - mV[,2])) 
  Cdiff_neg_samples <- as.numeric(Gamma_neg %*% sigma_mneg_vec)   # Samples from sqrt(n) (hat C_n - C)
  
  
  density_est_pos <- density(C+Cdiff_pos_samples/sqrt(nn), n=4096, from=-1, to=1)
  density_est_neg <- density(C+Cdiff_neg_samples/sqrt(nn), n=4096, from=-1, to=1)

  
  # ToDo: Add more densities!
  c_seq_long <- seq(-1, 1, length.out=10^4 + 1)
  df_density_hlp <- tibble(c_seq_long=c_seq_long,
                           i_set=i_set) %>%
    mutate(setting = paste0("p = ", p, ", q = ", q, ", C = ",C),
           setting_pq = paste0("p=", p, ", q=", q),
           setting_C = paste0("C=", C),
           type=type,
           density_norm_pos = dnorm(c_seq_long, mean=C, sd=sqrt(C_Var_pos)),
           density_norm_neg = dnorm(c_seq_long, mean=C, sd=sqrt(C_Var_neg)),
           density_halfnorm_zero = 1*(c_seq_long>=0) * dnorm(c_seq_long, mean=0, sd=sqrt(C_Var_zero_pos)) +
             1*(c_seq_long<0) * dnorm(c_seq_long, mean=0, sd=sqrt(C_Var_zero_neg)),
           density_equal_pos = approx(density_est_pos$x, density_est_pos$y, xout=c_seq_long)$y,
           density_equal_neg = approx(density_est_neg$x, density_est_neg$y, xout=c_seq_long)$y) %>%
    mutate(density = 
             dplyr::case_when(
               type == "normal_pos" ~ density_norm_pos,
               type == "normal_neg" ~ density_norm_neg,
               type == "zero" ~ density_halfnorm_zero,
               type == "equal_pos" ~ density_equal_pos,
               type == "equal_neg" ~ density_equal_neg
             ))
  
  df_density <- bind_rows(df_density, df_density_hlp)
}




# Generate the plot
df_plot <- res_df_MC %>%
  filter(measure=="C") %>%
  mutate(setting = paste0("p = ", p, ", q = ", q, ", C = ", Population),
         setting_pq = paste0("p=", p, ", q=", q),
         setting_C = paste0("C=", Population))

df_density_short <- full_join(df_density, 
          df_plot %>% 
            group_by(setting) %>%
            summarize(min=min(Estimate), max=max(Estimate)),
          by="setting") %>%
  dplyr::filter(c_seq_long >=min, c_seq_long<=max) %>%
  mutate(type=factor(type))


# Change factor levels
levels(df_density_short$type) <- list("Gaussian with C > 0    " = "normal_pos",        
                                      "Gaussian with C < 0    " = "normal_neg",
                                      "C = 0    " = "zero",
                                      "C > 0 and p = q    " = "equal_pos",
                                      "C < 0 and p = 1-q    " = "equal_neg")


ggplot(df_plot %>% filter(measure=="C")) +
  geom_histogram(aes(Estimate, after_stat(density)), boundary=0, bins=30, col="black", fill="gray") +
  geom_line(data=df_density_short, aes(x=c_seq_long, y=density, col=type), linewidth=1) +
  geom_vline(aes(xintercept=Population), linetype="dashed", linewidth=1, col="blue") +
  facet_wrap(~setting, scales="free") +
  xlab("Estimated C") +
  labs(col='Asymptotic distribution    ') +
  theme_bw() +
  theme(legend.position = "bottom") 


ggsave("simulation/output/C_AsyDist_new.pdf", width=24, height=16, unit="cm")

