## Test inference on simulated data
library(ViralLoadRJMCMC)
library(dplyr)

setwd("~/Documents/Research/Within-Host/RJMCMC_Results")

seeds <- 1:50

print(paste0("Analysis for seed", seeds[1]))

base_dir <- "~/Documents/Research/Within-Host/RJMCMC_Results/model_switch_results_good/"
# base_dir <- "~/Documents/Research/Within-Host/RJMCMC_Results/model_no_switch_results_unbounded_tp/"
# base_dir <- "~/Documents/Research/Within-Host/RJMCMC_Results/"

setwd(base_dir)

p_model_correct_all <- rep(NA,length(seeds))

mcmc_seed <- seeds[1]

set.seed(mcmc_seed)

analysis_dir <- paste0("./analysis_seed_",mcmc_seed)
output_dir <- paste0(analysis_dir,"/output")

setwd(analysis_dir) 

load("./quantile_coverage.RData")
load("./model_correct_summary.RData")
p_model_correct_all[1] <- p_model_correct

quantile_cover_df <- quantile_cover

for(seed_i in 2:length(seeds)){
  print(paste0("Analysis for seed", seeds[seed_i]))
  setwd(base_dir)
  
  mcmc_seed <- seeds[seed_i]
  
  set.seed(mcmc_seed)
  
  analysis_dir <- paste0("./analysis_seed_",mcmc_seed)
  output_dir <- paste0(analysis_dir,"/output")
 
  setwd(analysis_dir) 
  
  load("./quantile_coverage.RData")
  load("./model_correct_summary.RData")
  
  p_model_correct_all[seed_i] <- p_model_correct
  
  quantile_cover_df <- rbind(quantile_cover_df, quantile_cover)
}

quantile_cover_df$total_correct <- rowSums(quantile_cover_df)
quantile_cover_df$seed <- seeds
quantile_cover_df$p_model_correct <- p_model_correct_all