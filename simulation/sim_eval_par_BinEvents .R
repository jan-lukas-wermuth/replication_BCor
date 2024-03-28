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

# Load simulation data
res_df_MC <- readRDS(file = "simulation/data/sim_CI_20240308.rds")
alpha <- 0.9
M_min <- 500
n_plot <- c(100,500,2000)

pq_settings_filter <- c("p = 0.05,  q = 0.5", 
                        "p = 0.1,  q = 0.9", 
                        "p = 0.2,  q = 0.2", 
                        "p = 0.3,  q = 0.7", 
                        "p = 0.4,  q = 0.8") 


# Evaluation
eps <- 10^(-8)  # Use for numerical inaccuracies when comparing e.g. "0 <= 0"

# Compute coverage rates
CIcoverages <- res_df_MC %>% 
  # na.omit() %>%
  tidyr::drop_na(CI_lower,CI_upper) %>%
  dplyr::select(Estimate, CI_lower, CI_upper, i_MC, r, p, q, n, alpha, Fisher_Trans, measure, Population) %>%
  group_by(p,q,r,n,alpha,measure,Fisher_Trans) %>%
  summarize(coverage = mean(Population <= CI_upper + eps & Population >= CI_lower - eps),
            Lower = mean(Population <= CI_lower - eps),
            Upper = mean(Population >= CI_upper + eps),
            length = mean(CI_upper - CI_lower),
            pq_setting = paste0("p = ", first(p),",  q = ", first(q)),
            Population=round(first(Population),3),
            M=n()) %>%
  ungroup() %>%
  pivot_longer(cols=c(Lower, Upper), 
               names_to="CI Violation", values_to="miscoverage")


# Cosmetic changed for plotting
CIcoverages_plot <- CIcoverages %>%
  filter(M >= M_min,
         n %in% n_plot,
         pq_setting %in% pq_settings_filter) %>%
  rename(Asymptotics = Fisher_Trans) %>%
  rowwise() %>%
  mutate(Asymptotics = case_when(isTRUE(Asymptotics) ~ "Fisher Transformation",
                                 .default = "Standard"),
         Bound_lower = case_when(measure=="Phi" ~ (max(0,p+q-1)-p*q) / sqrt(p*(1-p)*q*(1-q)),
                                 .default = -1),
         Bound_upper = case_when(measure=="Phi" ~ (min(p,q)-p*q) / sqrt(p*(1-p)*q*(1-q)),
                                 .default = 1)) %>%
  arrange(desc(Asymptotics))



plot_height <- 12.5
plot_width <- 25


n_names <- c(`100` = "n = 100",
             `250` = "n = 250",
             `500` = "n = 500",
             `1000` = "n = 1000",
             `2000` = "n = 2000")


for (measure_plot in c("Phi", "C", "Q")){
  
  # Coverage rates w/o Fisher-Transformation
  p_coverage <- ggplot(CIcoverages_plot %>% 
         dplyr::filter(measure==measure_plot)) +
    geom_hline(yintercept=c(alpha), col="black") +
    geom_vline(aes(xintercept=Bound_lower), col="black", linetype="dotted") +
    geom_vline(aes(xintercept=Bound_upper), col="black", linetype="dotted") +
    geom_line(aes(x=Population, y=coverage, col=Asymptotics)) +
    scale_color_manual(breaks=c('Standard', 'Fisher Transformation'), values=c("red", "blue")) +
    coord_cartesian(ylim=c(0.7,1)) +
    facet_grid(n~pq_setting, labeller = labeller(n=n_names)) +
    ylab("Empirical Coverage Rate") + 
    xlab(paste0("True ", measure_plot)) +
    theme_bw() +
    theme(legend.position="bottom")

  
  # Directional (lower/upper) violation rates w/o Fisher-Transformation
  p_violation <- ggplot(CIcoverages_plot %>% 
           dplyr::filter(measure==measure_plot)) +
    geom_hline(yintercept=(1-alpha)/2, col="black") +
    geom_vline(aes(xintercept=Bound_lower), col="black", linetype="dotted") +
    geom_vline(aes(xintercept=Bound_upper), col="black", linetype="dotted") +
    geom_line(aes(x=Population, y=miscoverage, col=Asymptotics, linetype=`CI Violation`)) +
    scale_color_manual(breaks=c('Standard', 'Fisher Transformation'), values=c("red", "blue")) +
    coord_cartesian(ylim=c(0,0.2)) +
    facet_grid(n~pq_setting, labeller = labeller(n=n_names)) +
    ylab("Empirical Upper/Lower CI Violation Rate") + 
    xlab(paste0("True ", measure_plot)) +
    theme_bw() +
    theme(legend.position="bottom")
  
  
  # Delete legend for Q and C
  if (measure_plot %in% c("C", "Q")){
    p_coverage <- p_coverage + theme(legend.position="none")
    p_violation <- p_violation + theme(legend.position="none")
    
    plot_height_adj <- plot_height - 1.5
  } else {
    plot_height_adj <- plot_height
  }
  
  ggsave(paste0("simulation/output/", measure_plot, "_coverage.pdf"), plot=p_coverage,
         width=plot_width, height=plot_height_adj, unit="cm")
  
  ggsave(paste0("simulation/output/", measure_plot, "_violation.pdf"), plot=p_violation,
         width=plot_width, height=plot_height_adj, unit="cm")
}
