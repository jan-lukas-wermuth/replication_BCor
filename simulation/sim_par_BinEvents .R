# ===========================================================
# Possibly install the R package
# ===========================================================

# library(devtools)
# install_github("jan-lukas-wermuth/BCor")


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


# Specify the parameters
M <- 1000
n_set <- c(100, 500, 2000)
r_seq_length <- 41

alpha <- 0.9

# Choices for p and q
pq_list <- list(c(0.2, 0.2),
                c(0.3, 0.7),
                c(0.05, 0.5),
                c(0.4, 0.8),
                c(0.1, 0.1),
                c(0.1, 0.12),
                c(0.1, 0.9),
                c(0.05, 0.95),
                c(0.5, 0.5),
                c(0.01, 0.95))

# Loop in parallel
core.max <- 11
start.time <- Sys.time()
cl <- makeCluster(min(parallel::detectCores()-1, M, core.max) )
registerDoParallel(cl)
res_df_MC <- foreach(
  i_MC = 1:M,
  .combine=rbind,
  .packages=c("tidyverse", "tibble", "foreach", "BCor"),
  .errorhandling="pass"
)%dopar%{
  set.seed(i_MC) # set seed for reproducibility
  
  BGov_tbl <- tibble()
  tbl_coding <- tibble(X=c(1,0,1,0), Y=c(1,1,0,0))
  
  for (nn in n_set){
    for (i_pq in 1:length(pq_list)){
      p <- pq_list[[i_pq]][1]
      q <- pq_list[[i_pq]][2]
      
      # Equidistant "r_seq"
      r_min <- max(0,p+q-1)
      r_max <- min(p,q)
      r_seq <- seq(r_min, r_max, length.out = r_seq_length) %>% head(-1) %>% tail(-1)
      
      for (rr in r_seq){
        # Simulate
        prob <- c(rr, q-rr, p-rr, 1-p-q+rr) %>% pmax(0) %>% pmin(1)
        
        sample_index <- sample(1:4, size=nn, prob=prob,  replace = TRUE)
        sample_tbl <- tbl_coding[sample_index,]
        
        # Reduce the accuracy for the C-CI due to computational reasons
        c_seq_reduced <- seq(-1, 1, length.out = 401) %>% head(-1) %>% tail(-1)
        
        for (Fisher_logical in c(FALSE, TRUE)){
          
          Q_est <- YuleQ(sample_tbl$X, sample_tbl$Y, alpha=alpha, Fisher=Fisher_logical)
          Phi_est <- Phi(sample_tbl$X, sample_tbl$Y, alpha=alpha, Fisher=Fisher_logical)
          
          ########################################################################
          # Obtain details for the p-values for the estimated Cole coefficient
          C_est_details <- Cole(sample_tbl$X, sample_tbl$Y, alpha=alpha, Fisher=Fisher_logical,
                                m_rep = 10000, c_seq = c_seq_reduced, details_inference = TRUE)

          if (is_tibble(C_est_details)){ 
            C_est <- C_est_details %>% mutate(CI_lower=NA, CI_upper=NA) 
          } else {
            C_est <- C_est_details[[1]]
            C_est_pvals <- C_est_details[[2]] %>%
              mutate(pval_manual = pmin(1, 
                                        2*pmin(pmax(pv_uneq, 
                                                    2*pmin(pv_eq, 
                                                           pv_eq_pretest, 
                                                           na.rm = TRUE), 
                                                    na.rm = TRUE), 
                                               pv_pretest,
                                               na.rm = TRUE), 
                                        na.rm = TRUE),
                     pval_NoPQTest = pmin(1, 
                                          2*pmin(pmax(pv_uneq, 
                                                      pv_eq,
                                                      na.rm = TRUE), 
                                                 pv_pretest,
                                                 na.rm = TRUE),
                                          na.rm = TRUE),
                     pval_NoSigmaTest  = pmin(1, 
                                              pmax(pv_uneq, 
                                                   2*pmin(pv_eq, 
                                                          pv_eq_pretest, 
                                                          na.rm = TRUE), 
                                                   na.rm = TRUE), 
                                              na.rm = TRUE),
                     pval_NoPreTest = pmax(pv_uneq, 
                                           pv_eq,
                                           na.rm = TRUE))
            
            CI_manual <- C_est_pvals %>% dplyr::filter(pval_manual > (1-alpha)) %>% slice(c(1,n())) %>% pull(c_seq)
            CI_NoPQTest <- C_est_pvals %>% dplyr::filter(pval_NoPQTest > (1-alpha)) %>% slice(c(1,n())) %>% pull(c_seq)
            CI_NoSigmaTest <- C_est_pvals %>% dplyr::filter(pval_NoSigmaTest > (1-alpha)) %>% slice(c(1,n())) %>% pull(c_seq)
            CI_NoPreTest <- C_est_pvals %>% dplyr::filter(pval_NoPreTest > (1-alpha)) %>% slice(c(1,n())) %>% pull(c_seq)
            
            C_est <- as_tibble(C_est) %>% 
              mutate(CI_lower_manual = CI_manual[1],
                     CI_upper_manual = CI_manual[2],
                     CI_lower_NoPQTest = CI_NoPQTest[1],
                     CI_upper_NoPQTest = CI_NoPQTest[2],
                     CI_lower_NoSigmaTest = CI_NoSigmaTest[1],
                     CI_upper_NoSigmaTest = CI_NoSigmaTest[2],
                     CI_lower_NoPreTest = CI_NoPreTest[1],
                     CI_upper_NoPreTest = CI_NoPreTest[2])
          }
          ########################################################################
          
          if (!("CI_lower" %in% names(Q_est))){ Q_est <- Q_est %>% mutate(CI_lower=NA, CI_upper=NA) }
          if (!("CI_lower" %in% names(Phi_est))){ Phi_est <- Phi_est %>% mutate(CI_lower=NA, CI_upper=NA) }
          
          # Mutate tibbles to save them with all information
          C_est_mod <- as_tibble(C_est) %>% 
            rename(Estimate=C) %>% 
            mutate(i_MC=i_MC, r=rr, p=p, q=q, n=nn, alpha=alpha, Fisher_Trans=Fisher_logical,
                   measure="C", 
                   Population=1*(r-p*q >= 0) * (r-p*q)/(min(p,q) - p*q) + 1*(r-p*q < 0) * (r-p*q)/(p*q - max(0, p+q-1)))
          
          Q_est_mod <- as_tibble(Q_est) %>% 
            rename(Estimate=Q) %>% 
            mutate(i_MC=i_MC, r=rr, p=p, q=q, n=nn, alpha=alpha, Fisher_Trans=Fisher_logical,
                   measure="Q", 
                   Population=(r-p*q) / (r * (1-p-q+r) + (q-r) * (p-r)))
          
          Phi_est_mod <- as_tibble(Phi_est) %>% 
            rename(Estimate=Phi) %>% 
            mutate(i_MC=i_MC, r=rr, p=p, q=q, n=nn, alpha=alpha, Fisher_Trans=Fisher_logical,
                   measure="Phi", 
                   Population=(r-p*q) / (sqrt(p*(1-p)*q*(1-q))))
          
          
          # Merge
          BGov_tbl <- bind_rows(BGov_tbl,
                                bind_rows(C_est_mod, Q_est_mod, Phi_est_mod))
          
        }
      }
    }
  }
  
  BGov_tbl
}  
stopCluster(cl)
end.time <- Sys.time()
(run.time <- end.time-start.time)

head(res_df_MC)
saveRDS(res_df_MC, file = "simulation/data/sim_CI_20240308.rds")


