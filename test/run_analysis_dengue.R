## Run an RJMCMC analysis
# Read in clean data
library(dplyr)

library(ViralLoadRJMCMC)

setwd("~/Documents/Research/Within-Host/RJMCMC_Results")

mcmc_seed <- 1111
set.seed(mcmc_seed)

analysis_dir <- paste0("~/Documents/Research/Within-Host/RJMCMC_Results/analysis_dengue_seed_",mcmc_seed)
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
load("~/Documents/Research/Within-Host/Within_Host_CT_Analysis/Data/Dengue/clean_dengue.RData")

clean_data$index <- match(clean_data$No,sort(unique(clean_data$No))) - 1

indiv_dataset <- clean_data %>%
  arrange(index) %>%
  group_by(No) %>%
  filter(l_value == max(l_value)) %>%
  mutate(dp_min = l_value - 4) %>%
  mutate(max_viral_load = l_value) %>%
  select(No,index,plasma_leakage,plasma_leakage_n,dp_min,max_viral_load,first_gt0_adj,last_gt0_adj,min_obs_day_adj,max_obs_day_adj,n_gt0,zero_begin,zero_end,est_wp,est_wp_wr,est_wr)

init_model <- ifelse(indiv_dataset$est_wp,1,
                     ifelse(indiv_dataset$est_wp_wr,2,
                            ifelse(indiv_dataset$est_wr,3,0)))

infect_data <- clean_data

## Change subtype to 0 and 1 and 2
n_subtype <- length(unique(infect_data$plasma_leakage_n))
infect_data$subtype_inf_n_index <- (0:(n_subtype-1))[infect_data$plasma_leakage_n]
indiv_dataset$subtype_inf_n_index <- (0:(n_subtype-1))[indiv_dataset$plasma_leakage_n]

infect_data <- infect_data %>%
                arrange(index,day_adj)

viral_data <- data.frame(index=infect_data$index,
                         viral_load=infect_data$l_value,
                         time=infect_data$day_adj)

individual_data <- data.frame(subtype=indiv_dataset$subtype_inf_n_index,
                              t_first_positive=indiv_dataset$first_gt0_adj,
                              t_last_positive=indiv_dataset$last_gt0_adj,
                              t_first_test=indiv_dataset$min_obs_day_adj,
                              t_last_test=indiv_dataset$max_obs_day_adj,
                              n_positive=indiv_dataset$n_gt0,
                              max_viral_load=indiv_dataset$max_viral_load)

settings <- data.frame(lod=20,
                       sensitivity=0.99,
                       n_iterations=200000,
                       n_subtypes=n_subtype,
                       n_subjects=nrow(individual_data),
                       n_data=nrow(viral_data))

priors <- data.frame(
  wp_min = 0.5,
  wp_max = 21,
  wr_min = 0.5,
  wr_max = 25,
  wpmean_max = 14,
  dpmean_max = 40,
  wrmean_max = 25,
  wpsd_max = 5,
  tpsd_max = 5,
  dpsd_max = 10,
  wrsd_max = 5,
  sigma_max = 15,
  wpsd_min = 0.3,
  tpsd_min = 0.3,
  dpsd_min = 0.3,
  wrsd_min = 0.3,
  sigma_min = 0,
  wpmean_mean = 5,
  wpmean_sd = 100,
  dpmean_mean = 15,
  dpmean_sd = 100,
  wrmean_mean = 6,
  wrmean_sd = 100,
  wpsd_scale = 100,
  tpsd_scale = 100,
  dpsd_scale = 100,
  wrsd_scale = 100,
  sigma_scale = 100
)

wp_mean_init <- rep(2,settings$n_subtypes)
wp_sd_init <- rep(1,settings$n_subtypes)
dp_mean_init <- rep(15,settings$n_subtypes)
dp_sd_init <- rep(5,settings$n_subtypes)
tp_sd_init <- rep(2,settings$n_subtypes)
wr_mean_init <- rep(5,settings$n_subtypes)
wr_sd_init <- rep(1,settings$n_subtypes)
wp_init <- runif(settings$n_subjects,1,9)
dp_init <- rnorm(settings$n_subjects,15,3)
tp_init <- rnorm(settings$n_subjects,0,1)
wr_init <- runif(settings$n_subjects,1,8)
model_init <- ifelse(indiv_dataset$est_wp,1,
                     ifelse(indiv_dataset$est_wp_wr,2,
                            ifelse(indiv_dataset$est_wr,3,0)))

sigma_init <- 5

wp_mean_sf <- rep(1,settings$n_subtypes)
wp_sd_sf <- rep(0.5,settings$n_subtypes)
dp_mean_sf <- rep(4,settings$n_subtypes)
dp_sd_sf <- rep(2,settings$n_subtypes)
tp_sd_sf <- rep(0.5,settings$n_subtypes)
wr_mean_sf <- rep(1,settings$n_subtypes)
wr_sd_sf <- rep(0.5,settings$n_subtypes)
sigma_sf <- 0.2
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

source("~/Documents/Research/Within-Host/ViralLoadRJMCMC/test/summarize_results_dengue.R")