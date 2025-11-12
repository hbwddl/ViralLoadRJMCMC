## Run an RJMCMC analysis
# Read in clean data
library(dplyr)

library(ViralLoadRJMCMC)

setwd("~/Documents/Research/Within-Host/RJMCMC_Results")

mcmc_seed <- 1115
set.seed(mcmc_seed)

analysis_dir <- paste0("~/Documents/Research/Within-Host/RJMCMC_Results/analysis_flu_seed_",mcmc_seed)
output_dir <- paste0(analysis_dir,"/output")

if(!dir.exists(analysis_dir)){
  dir.create(analysis_dir)
  setwd(analysis_dir)
} else{
  setwd(analysis_dir)
}

if(!dir.exists(output_dir)){
  dir.create(output_dir)
  setwd(output_dir)
} else{
  setwd(output_dir)
}

setwd(analysis_dir)

## Get data
# Load clean data and individual data
load("~/Documents/Research/Within-Host/Within_Host_CT_Analysis/Data/clean_data.RData")
load("~/Documents/Research/Within-Host/Within_Host_CT_Analysis/Data/indiv_dataset.RData")

infect_data <- clean_data %>%
                filter(subtype_inf_n != 0) %>%
                filter(pig_id %in% indiv_dataset$pig_id)

## Change subtype to 0 and 1 and 2
n_subtype <- length(unique(infect_data$subtype_inf_n))
infect_data$subtype_inf_n_index <- (0:(n_subtype-1))[infect_data$subtype_inf_n]
indiv_dataset$subtype_inf_n_index <- (0:(n_subtype-1))[indiv_dataset$subtype_inf_n]

infect_data$index <- indiv_dataset$index[match(infect_data$pig_id,indiv_dataset$pig_id)]

infect_data <- infect_data %>%
                arrange(index,day_adj)

infect_data_filter <- infect_data %>%
                      filter(index != 147)

indiv_dataset_filter <- indiv_dataset %>%
                        filter(index != 147)

indiv_dataset_filter$index_adj <- 0:(nrow(indiv_dataset_filter)-1)

infect_data_filter$index_adj <- indiv_dataset_filter$index_adj[match(infect_data_filter$index,indiv_dataset_filter$index)]

viral_data <- data.frame(index=infect_data_filter$index_adj,
                         viral_load=infect_data_filter$ct_delta,
                         time=infect_data_filter$day_adj)

individual_data <- data.frame(index=indiv_dataset_filter$index_adj,
                              subtype=indiv_dataset_filter$subtype_inf_n_index,
                              t_first_positive=indiv_dataset_filter$first_gt0_adj,
                              t_last_positive=indiv_dataset_filter$last_gt0_adj,
                              t_first_test=indiv_dataset_filter$min_obs_day_adj,
                              t_last_test=indiv_dataset_filter$max_obs_day_adj,
                              n_positive=indiv_dataset_filter$n_gt0,
                              max_viral_load=indiv_dataset_filter$max_viral_load)

settings <- data.frame(lod=45,
                       sensitivity=1,
                       n_iterations=500000,
                       n_subtypes=n_subtype,
                       n_subjects=nrow(individual_data),
                       n_data=nrow(viral_data))

priors <- data.frame(
  wp_min = 0.5,
  wp_max = 10,
  wr_min = 0.5,
  wr_max = 25,
  wpmean_max = 8,
  dpmean_max = 40,
  wrmean_max = 25,
  wpsd_max = 10,
  tpsd_max = 10,
  dpsd_max = 10,
  wrsd_max = 10,
  sigma_max = 20,
  wpsd_min = 0.3,
  tpsd_min = 0.3,
  dpsd_min = 0.3,
  wrsd_min = 0.3,
  sigma_min = 0,
  wpmean_mean = 5,
  wpmean_sd = 3,
  dpmean_mean = 30,
  dpmean_sd = 100,
  wrmean_mean = 7,
  wrmean_sd = 3,
  wpsd_scale = 3,
  tpsd_scale = 100,
  dpsd_scale = 100,
  wrsd_scale = 3,
  sigma_scale = 100
)

wp_mean_init <- rep(6,settings$n_subtypes)
wp_sd_init <- rep(1,settings$n_subtypes)
dp_mean_init <- rep(25,settings$n_subtypes)
dp_sd_init <- rep(10,settings$n_subtypes)
tp_sd_init <- rep(2,settings$n_subtypes)
wr_mean_init <- rep(6,settings$n_subtypes)
wr_sd_init <- rep(1,settings$n_subtypes)
wp_init <- rep(5,settings$n_subjects)
dp_init <- rep(25,settings$n_subjects)
tp_init <- rnorm(settings$n_subjects,0,1)
tp_init <- rep(0,settings$n_subjects)
wr_init <- rep(6,settings$n_subjects)
model_init <- sample(c(1,2,3),settings$n_subjects,replace = T)
model_init <- rep(2,settings$n_subjects)
model_init[1:4] <- 1
model_init[9:16] <- 1
model_init[23] <- 1
model_init[34] <- 1
model_init[38:42] <- 1
model_init[48:50] <- 1
model_init[75:78] <- 1
model_init[82:83] <- 1
model_init[90:92] <- 1
model_init[94] <- 1
model_init[144] <- 1
model_init[148] <- 3
model_init[68] <- 1
model_init[133] <- 1
sigma_init <- 5

wp_mean_sf <- rep(2,settings$n_subtypes)
wp_sd_sf <- rep(1.75,settings$n_subtypes)
dp_mean_sf <- rep(5,settings$n_subtypes)
dp_sd_sf <- rep(4,settings$n_subtypes)
tp_sd_sf <- rep(0.5,settings$n_subtypes)
wr_mean_sf <- rep(1.5,settings$n_subtypes)
wr_sd_sf <- rep(1.5,settings$n_subtypes)
sigma_sf <- 0.1
wp_sf <- 1
tp_sf <- 1
dp_sf <- 1
wr_sf <- 1

save(viral_data,file="viral_data_in.RData")
save(individual_data,file="individual_data_in.RData")
save(settings,mcmc_seed,file="settings_in.RData")
save(priors,file="priors.RData")

sink(file="rjmcmc.out")

viral_load_rjmcmc("./output/",
                  viral_data,
                  individual_data,
                  settings,
                  wp_mean_init,
                  wp_sd_init,
                  dp_mean_init,
                  dp_sd_init,
                  tp_sd_init,
                  wr_mean_init,
                  wr_sd_init,
                  wp_init,
                  dp_init,
                  tp_init,
                  wr_init,
                  model_init,
                  sigma_init,
                  priors,
                  wp_mean_sf,
                  wp_sd_sf,
                  dp_mean_sf,
                  dp_sd_sf,
                  tp_sd_sf,
                  wr_mean_sf,
                  wr_sd_sf,
                  sigma_sf,
                  wp_sf,
                  tp_sf,
                  dp_sf,
                  wr_sf,
                  mcmc_seed)

sink(file=NULL)

source("~/Documents/Research/Within-Host/ViralLoadRJMCMC/test/summarize_results_flu.R")