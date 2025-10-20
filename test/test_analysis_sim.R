## Test inference on simulated data
library(ViralLoadRJMCMC)
library(dplyr)

setwd("~/Documents/Research/Within-Host/RJMCMC_Results")

mcmc_seed <- 1212
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
}

data_settings_in <- list(lod=45,
                         n=400,
                         p_group=c(0.33,0.34,0.33),
                         t_obs=1:7,
                         sensitivity=1)

param_settings_in <- list(p_model=c(0,0.6,0.4),
                          wp_mean=c(4,4,4),
                          wp_sd=c(1,1,1),
                          tp_sd=c(1,1,1),
                          dp_mean=c(30,30,30),
                          dp_sd=c(5,5,5),
                          wr_mean=c(6,6,6),
                          wr_sd=c(1,1,1),
                          sigma=5,
                          wp_min=1,
                          wp_max=8,
                          dp_min=20,
                          tp_min=-2,
                          tp_max=2,
                          wr_min=1,
                          wr_max=10)

save(param_settings_in,file="param_settings_in.RData")

sim_out <- simulate_viral_load_data(data_settings_in,
                                    param_settings_in,
                                    sim_seed=mcmc_seed)

plot_sim_data(sim_out)

save(sim_out,file="sim_out.RData")

n_subtype <- length(sim_out$settings$p_group)

individual_data <- data.frame(index_init=sim_out$individual_data$index,
                              subtype=sim_out$individual_data$group,
                              t_first_positive=sim_out$individual_data$t_first_positive,
                              t_last_positive=sim_out$individual_data$t_last_positive,
                              t_first_test=sim_out$individual_data$t_first_test,
                              t_last_test=sim_out$individual_data$t_last_test,
                              n_positive=sim_out$individual_data$n_positive,
                              max_viral_load=sim_out$individual_data$max_viral_load,
                              model_true=sim_out$individual_data$model,
                              tp_true=sim_out$individual_data$tp,
                              wp_true=sim_out$individual_data$wp,
                              dp_true=sim_out$individual_data$dp,
                              wr_true=sim_out$individual_data$wr) %>%
                              filter(n_positive > 1)

viral_data <- data.frame(index_init=sim_out$viral_load_data$index,
                         viral_load=sim_out$viral_load_data$viral_load_obs,
                         time=sim_out$viral_load_data$time) %>%
                        filter(index_init %in% individual_data$index_init)

individual_data$index <- 0:(nrow(individual_data)-1)
viral_data$index <- match(viral_data$index_init,individual_data$index_init)-1

settings <- data.frame(lod=45,
                       sensitivity=0.99,
                       n_iterations=100000,
                       n_subtypes=n_subtype,
                       n_subjects=nrow(individual_data),
                       n_data=nrow(viral_data))

priors <- data.frame(wp_min = 0.5,
                      wp_max = 20,
                      wr_min = 0.5,
                      wr_max = 20,
                      wpmean_max = 20,
                      dpmean_max = 40,
                      wrmean_max = 20,
                      wpsd_max = 10,
                      tpsd_max = 10,
                      dpsd_max = 20,
                      wrsd_max = 10,
                      sigma_max = 20,
                      wpsd_min = 0,
                      tpsd_min = 0,
                      dpsd_min = 0,
                      wrsd_min = 0,
                      sigma_min = 0,
                      wpmean_mean = 5,
                      wpmean_sd = 100,
                      dpmean_mean = 30,
                      dpmean_sd = 100,
                      wrmean_mean = 7,
                      wrmean_sd = 5,
                      wpsd_scale = 100,
                      tpsd_scale = 100,
                      dpsd_scale = 100,
                      wrsd_scale = 100,
                      sigma_scale = 100
)

wp_mean_init <- rep(3,settings$n_subtypes)
wp_sd_init <- rep(1,settings$n_subtypes)
dp_mean_init <- rep(25,settings$n_subtypes)
dp_sd_init <- rep(10,settings$n_subtypes)
tp_sd_init <- rep(2,settings$n_subtypes)
wr_mean_init <- rep(5,settings$n_subtypes)
wr_sd_init <- rep(1,settings$n_subtypes)
wp_init <- runif(settings$n_subjects,1,9)
dp_init <- rnorm(settings$n_subjects,25,10)
tp_init <- rnorm(settings$n_subjects,0,2)
wr_init <- runif(settings$n_subjects,1,8)
model_init <- rep(2,settings$n_subjects)
sigma_init <- 10

wp_mean_init <- sim_out$parameters$wp_mean
wp_sd_init <- sim_out$parameters$wp_sd
dp_mean_init <- sim_out$parameters$dp_mean
dp_sd_init <- sim_out$parameters$dp_sd
tp_sd_init <- sim_out$parameters$tp_sd
wr_mean_init <- sim_out$parameters$wr_mean
wr_sd_init <- sim_out$parameters$wr_sd
wp_init <- individual_data$wp_true
dp_init <- individual_data$dp_true
tp_init <- individual_data$tp_true
wr_init <- individual_data$wr_true
model_init <- individual_data$model_true
sigma_init <- sim_out$parameters$sigma

wp_mean_sf <- rep(2,settings$n_subtypes)
wp_sd_sf <- rep(0.5,settings$n_subtypes)
dp_mean_sf <- rep(5,settings$n_subtypes)
dp_sd_sf <- rep(4,settings$n_subtypes)
tp_sd_sf <- rep(0.5,settings$n_subtypes)
wr_mean_sf <- rep(1.5,settings$n_subtypes)
wr_sd_sf <- rep(0.5,settings$n_subtypes)
sigma_sf <- 0.1
wp_sf <- 0.5
tp_sf <- 0.5
dp_sf <- 1
wr_sf <- 0.5

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

source("~/Documents/Research/Within-Host/ViralLoadRJMCMC/test/summarize_results_sim.R")