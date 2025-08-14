## Run an RJMCMC analysis
# Read in clean data
library(dplyr)

library(ViralLoadRJMCMC)

setwd("~/Documents/Research/Within-Host/RJMCMC_Results")

mcmc_seed <- 1111
set.seed(mcmc_seed)

analysis_dir <- paste0("~/Documents/Research/Within-Host/RJMCMC_Results/analysis_seed_",mcmc_seed)
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

viral_data <- data.frame(index=infect_data$index,
                         viral_load=infect_data$ct_delta,
                         time=infect_data$day_adj)

individual_data <- data.frame(subtype=indiv_dataset$subtype_inf_n_index,
                              t_first_positive=indiv_dataset$first_gt0_adj,
                              t_last_positive=indiv_dataset$last_gt0_adj,
                              t_first_test=indiv_dataset$min_obs_day_adj,
                              t_last_test=indiv_dataset$max_obs_day_adj,
                              n_positive=indiv_dataset$n_gt0)

settings <- data.frame(lod=45,
                       sensitivity=0.99,
                       n_iterations=1000,
                       n_subtypes=3,
                       n_subjects=203,
                       n_data=1092)

priors <- data.frame(
  wp_max = 20,
  wr_max = 20,
  wpmean_max = 20,
  dpmean_max = 40,
  wrmean_max = 20,
  wpsd_max = 10,
  tpsd_max = 10,
  dpsd_max = 20,
  wrsd_max = 10,
  sigma_max = 10,
  wpsd_min = 0,
  tpsd_min = 0,
  dpsd_min = 0,
  wrsd_min = 0,
  sigma_min = 0,
  wpmean_mean = 5,
  wpmean_sd = 100,
  dpmean_mean = 30,
  dpmean_sd = 100,
  wrmean_mean = 5,
  wrmean_sd = 100,
  wpsd_scale = 100,
  tpsd_scale = 100,
  dpsd_scale = 100,
  wrsd_scale = 100,
  sigma_scale = 100
)



viral_load_rjmcmc("./output/",
                  viral_data,
                  individual_data,
                  settings,
                  wp_mean_init_vec,
                  wp_sd_init_vec,
                  dp_mean_init_vec,
                  dp_sd_init_vec,
                  tp_sd_init_vec,
                  wr_mean_init_vec,
                  wr_sd_init_vec,
                  wp_init_vec,
                  dp_init_vec,
                  tp_init_vec,
                  wr_init_vec,
                  model_init_vec,
                  sigma_init,
                  priors,
                  wp_mean_sf_vec,
                  wp_sd_sf_vec,
                  dp_mean_sf_vec,
                  dp_sd_sf_vec,
                  tp_sd_sf_vec,
                  wr_mean_sf_vec,
                  wr_sd_sf_vec,
                  sigma_sf,
                  wp_sf,
                  tp_sf,
                  dp_sf,
                  wr_sf,
                  seed_arg)
