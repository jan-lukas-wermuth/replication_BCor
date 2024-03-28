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
library(patchwork)


# Load simulation data
res_df_MC <- readRDS(file = "simulation/data/sim_CI_20240308.rds")
alpha <- 0.9
M_min <- 500
n_plot <- c(500)

pq_settings_filter <- c("p = 0.05,  q = 0.5", 
                        "p = 0.1,  q = 0.9", 
                        "p = 0.2,  q = 0.2", 
                        "p = 0.3,  q = 0.7", 
                        "p = 0.4,  q = 0.8") 

# Evaluation
eps <- 10^(-8) # Use for numerical inaccuracies when comparing e.g. "0 <= 0"

res_df_MC <- res_df_MC %>%
  dplyr::filter(measure=="C") %>%
  rename(CI_lower_standard = CI_lower,
         CI_upper_standard = CI_upper)


res_df_MC_long <- full_join(res_df_MC %>% 
                              dplyr::select(!c(CI_upper_standard, CI_upper_manual, CI_upper_NoPQTest, CI_upper_NoSigmaTest, CI_upper_NoPreTest)) %>%    
                              pivot_longer(cols = starts_with("CI_lower"),
                                           names_to = "type",
                                           names_prefix = "CI_lower_",
                                           values_to = "lower"),
                            res_df_MC %>% 
                              dplyr::select(!c(CI_lower_standard, CI_lower_manual, CI_lower_NoPQTest, CI_lower_NoSigmaTest, CI_lower_NoPreTest)) %>%     
                              pivot_longer(cols = starts_with("CI_upper"),
                                           names_to = "type",
                                           names_prefix = "CI_upper_",
                                           values_to = "upper"),
                            by=c("Estimate", "i_MC", "r", "p", "q", "n", "alpha", "Fisher_Trans", "measure", "Population", "type"))
  
               
# Compute coverage rates
CIcoverages <- res_df_MC_long %>% 
  na.omit() %>%
  group_by(p,q,r,n,alpha,measure, Fisher_Trans, type) %>%
  summarize(coverage = mean(Population <= upper + eps & Population >= lower - eps),
            length = mean(upper - lower),
            Lower = mean(Population <= lower - eps),
            Upper = mean(Population >= upper + eps),
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
         pq_setting %in% pq_settings_filter
         ) %>%
  rename(Asymptotics = Fisher_Trans) %>%
  rowwise() %>%
  mutate(Asymptotics = case_when(isTRUE(Asymptotics) ~ "Fisher Transformation",
                                 .default = "Standard"),
         Bound_lower = case_when(measure=="Phi" ~ (max(0,p+q-1)-p*q) / sqrt(p*(1-p)*q*(1-q)),
                                 .default = -1),
         Bound_upper = case_when(measure=="Phi" ~ (min(p,q)-p*q) / sqrt(p*(1-p)*q*(1-q)),
                                 .default = 1)) %>%
  arrange(desc(Asymptotics))


# Filter some settings we want to plot
CIcoverages_plot2 <- CIcoverages_plot %>% 
  dplyr::filter(measure=="C",
                Asymptotics=="Fisher Transformation",
                type != "manual")


# Cosmetic changes in the plots
CIcoverages_plot2$type <- factor(CIcoverages_plot2$type, 
                                 levels=c("standard", "NoSigmaTest", "NoPQTest", "NoPreTest"), 
                                 labels = c("Full Combination", "No Sigma Test", "No p,q Test", "Basic Combination"))

colors_manual <- c("Full Combination"="red",
                   "No Sigma Test"="deepskyblue2",
                   "No p,q Test"="darkblue",
                   "Basic Combination"="black")

linetypes_manual <- c("Full Combination"="solid",
                      "No Sigma Test"="dashed",
                      "No p,q Test"="dashed",
                      "Basic Combination"="dotted")

# Coverage rates of different CI construction methods
p_coverage <- ggplot(CIcoverages_plot2) +
  geom_hline(yintercept=c(alpha), col="black") +
  geom_line(aes(x=Population, y=coverage, col=type, linetype=type)) +
  scale_color_manual(values=colors_manual, name="Confidence Interval Construction") +
  scale_linetype_manual(values=linetypes_manual, name="Confidence Interval Construction") +
  coord_cartesian(ylim=c(0.7,1)) +
  facet_grid(~pq_setting) +
  ylab("Empirical Coverage Rate") + 
  xlab(paste0("True C")) +
  theme_bw() +
  theme(legend.position="bottom")


# CI lengths of different CI construction methods
p_CIlength <- ggplot(CIcoverages_plot2) +
  geom_line(aes(x=Population, y=length, col=type, linetype=type)) +
  scale_color_manual(values=colors_manual, name="Confidence Interval Construction") +
  scale_linetype_manual(values=linetypes_manual, name="Confidence Interval Construction") +
  facet_grid(~pq_setting) +
  ylab("Average CI Length") + 
  xlab(paste0("True C")) +
  theme_bw() +
  theme(legend.position="bottom")

p_CIlength


p_C <- ( p_coverage / p_CIlength ) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
p_C

ggsave(paste0("simulation/output/C_details.pdf"), 
       plot=p_C,
       width=25, height=13.5, unit="cm")







