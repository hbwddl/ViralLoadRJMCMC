## Fit data to covid
library(dplyr)
library(ViralLoadRJMCMC)

setwd("~/Documents/Research/Within-Host/RJMCMC_Results")

load("~/Documents/Research/Within-Host/Within_Host_CT_Analysis/Data/ct_dat_analysis.RData")


mcmc_seed <- 1113
set.seed(mcmc_seed)

analysis_dir <- paste0("~/Documents/Research/Within-Host/RJMCMC_Results/analysis_covid_seed_",mcmc_seed)
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

analysis_data$index <- analysis_data$id_clean - 1

analysis_data_filter <- analysis_data %>%
                        # filter(!(id_clean %in% c(5,23,36))) %>%
                        group_by(id_clean) %>%
                        mutate(max_viral_load = max(ct_delta)) %>%
                        ungroup() 
# %>%
#                         filter(max_viral_load > 10)
                          
individual_data_select <- analysis_data_filter %>%
                    group_by(index) %>%
                    arrange(day_adj) %>%
                    mutate(n_gt0 = sum(is_gt0)) %>% 
                    mutate(max_viral_load = max(ct_delta)) %>%
                    slice(1) %>%
                    ungroup() %>%
                    # filter(n_gt0 > 1) %>%
                    arrange(index)

individual_data_select$index_adj <- 0:(nrow(individual_data_select)-1)

analysis_data_select <- analysis_data_filter %>%
                        filter(index %in% individual_data_select$index) %>%
                        arrange(index) %>%
                        filter(!(index == 18 & day_adj %in% c(2,4))) %>%
                        filter(!(index == 43 & day_adj %in% c(1,2,3))) %>%
                        # filter(!(index == 28 & day_adj %in% c(-8,-7,-2,-1))) %>%
                        filter(!(index == 4 & day_adj %in% c(11,12))) %>%
                        filter(!(index == 32 & day_adj %in% c(1,2)))

analysis_data_select$index_adj <- individual_data_select$index_adj[match(analysis_data_select$index,individual_data_select$index)]

n_subtype <- length(unique(analysis_data$symptomatic))

viral_data <- data.frame(index=analysis_data_select$index_adj,
                         viral_load=analysis_data_select$ct_delta,
                         time=analysis_data_select$day_adj)

individual_data <- data.frame(index=individual_data_select$index_adj,
                              subtype=individual_data_select$symptomatic,
                              t_first_positive=individual_data_select$first_gt0_adj,
                              t_last_positive=individual_data_select$last_gt0_adj,
                              t_first_test=individual_data_select$min_obs_day_adj,
                              t_last_test=individual_data_select$max_obs_day_adj,
                              n_positive=individual_data_select$n_gt0,
                              max_viral_load=individual_data_select$max_viral_load)

settings <- data.frame(lod=40,
                       sensitivity=0.99,
                       n_iterations=200000,
                       n_subtypes=n_subtype,
                       n_subjects=nrow(individual_data),
                       n_data=nrow(viral_data))

priors <- data.frame(
  wp_min = 0.5,
  wp_max = 20,
  wr_min = 0.5,
  wr_max = 25,
  wpmean_max = 10,
  dpmean_max = 40,
  wrmean_max = 20,
  wpsd_max = 10,
  tpsd_max = 10,
  dpsd_max = 10,
  wrsd_max = 10,
  sigma_max = 20,
  wpsd_min = 0.5,
  tpsd_min = 0.3,
  dpsd_min = 2,
  wrsd_min = 1,
  sigma_min = 1,
  wpmean_mean = 3,
  wpmean_sd = 100,
  dpmean_mean = 25,
  dpmean_sd = 100,
  wrmean_mean = 14,
  wrmean_sd = 100,
  wpsd_scale = 100,
  tpsd_scale = 100,
  dpsd_scale = 100,
  wrsd_scale = 100,
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
model_init <- rep(2,settings$n_subjects)
sigma_init <- 5

wp_mean_sf <- rep(2,settings$n_subtypes)
wp_sd_sf <- rep(1.75,settings$n_subtypes)
dp_mean_sf <- rep(5,settings$n_subtypes)
dp_sd_sf <- rep(4,settings$n_subtypes)
tp_sd_sf <- rep(0.5,settings$n_subtypes)
wr_mean_sf <- rep(1.5,settings$n_subtypes)
wr_sd_sf <- rep(1.5,settings$n_subtypes)
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

source("~/Documents/Research/Within-Host/ViralLoadRJMCMC/test/summarize_results_covid.R")