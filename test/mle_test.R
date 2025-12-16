## Test inference on simulated data
library(ViralLoadRJMCMC)
library(dplyr)

setwd("~/Documents/Research/Within-Host/RJMCMC_Results")

mcmc_seed <- 6
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
                         t_obs=(-10):10,
                         sensitivity=1)

param_settings_in <- list(p_model=c(0.45,0.45,0.1),
                          wp_mean=c(3.5,3.5,3.5),
                          wp_sd=c(1,1,1),
                          tp_sd=c(2,2,2),
                          dp_mean=c(24,24,24),
                          dp_sd=c(3,3,3),
                          wr_mean=c(7,12,9),
                          wr_sd=c(1,1,1),
                          sigma=4.5,
                          wp_min=0.5,
                          wp_max=20,
                          dp_min=15,
                          tp_min=-2,
                          tp_max=2,
                          wr_min=0.5,
                          wr_max=20)

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
                     sigma_min = 4,
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

save(viral_data,file="viral_data_in.RData")
save(individual_data,file="individual_data_in.RData")
save(settings,mcmc_seed,file="settings_in.RData")
save(priors,file="priors.RData")

log_likelihood_r(viral_data,
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
                  sigma_init)

sigma_val <- seq(from=0.5,to=10,length.out=200)
lh_val <- rep(NA,length(sigma_val))

for(i in 1:length(sigma_val)){
  lh_val[i] <- log_likelihood_r(viral_data,
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
                                sigma_val[i])
}

plot(sigma_val,lh_val,type="l",main="Sigma profile likelihood")
abline(v=sigma_init,col="blue")

print(sigma_val[which.max(lh_val)])

wr_mean_1_val <- seq(from=wr_mean_init[1]/2,to=wr_mean_init[1]*2,length.out=200)
wr_lh_1_val <- wr_mean_1_val
wr_mean_vec <- wr_mean_init

for(i in 1:length(sigma_val)){
  wr_mean_vec[1] <- wr_mean_1_val[i]
  wr_lh_1_val[i] <- log_likelihood_r(viral_data,
                                individual_data,
                                settings,
                                wp_mean_init,
                                wp_sd_init,
                                dp_mean_init,
                                dp_sd_init,
                                tp_sd_init,
                                wr_mean_vec,
                                wr_sd_init,
                                wp_init,
                                dp_init,
                                tp_init,
                                wr_init,
                                model_init,
                                sigma_init)
}

plot(wr_mean_1_val,wr_lh_1_val,type="l",main="WR Mean 1 profile likelihood")
abline(v=wr_mean_init[1],col="blue")
print(wr_mean_1_val[which.max(wr_lh_1_val)])

wr_mean_2_val <- seq(from=wr_mean_init[2]/2,to=wr_mean_init[2]*2,length.out=200)
wr_lh_2_val <- wr_mean_2_val
wr_mean_vec <- wr_mean_init

for(i in 1:length(sigma_val)){
  wr_mean_vec[2] <- wr_mean_2_val[i]
  wr_lh_2_val[i] <- log_likelihood_r(viral_data,
                                     individual_data,
                                     settings,
                                     wp_mean_init,
                                     wp_sd_init,
                                     dp_mean_init,
                                     dp_sd_init,
                                     tp_sd_init,
                                     wr_mean_vec,
                                     wr_sd_init,
                                     wp_init,
                                     dp_init,
                                     tp_init,
                                     wr_init,
                                     model_init,
                                     sigma_init)
}

plot(wr_mean_2_val,wr_lh_2_val,type="l",main="WR Mean 2 profile likelihood")
abline(v=wr_mean_init[2],col="blue")
print(wr_mean_2_val[which.max(wr_lh_2_val)])

wr_mean_3_val <- seq(from=wr_mean_init[3]/2,to=wr_mean_init[3]*2,length.out=200)
wr_lh_3_val <- wr_mean_3_val
wr_mean_vec <- wr_mean_init

for(i in 1:length(sigma_val)){
  wr_mean_vec[3] <- wr_mean_3_val[i]
  wr_lh_3_val[i] <- log_likelihood_r(viral_data,
                                     individual_data,
                                     settings,
                                     wp_mean_init,
                                     wp_sd_init,
                                     dp_mean_init,
                                     dp_sd_init,
                                     tp_sd_init,
                                     wr_mean_vec,
                                     wr_sd_init,
                                     wp_init,
                                     dp_init,
                                     tp_init,
                                     wr_init,
                                     model_init,
                                     sigma_init)
}

plot(wr_mean_3_val,wr_lh_3_val,type="l",main="WR Mean 3 profile likelihood")
abline(v=wr_mean_init[3],col="blue")
print(wr_mean_3_val[which.max(wr_lh_3_val)])

wr_sd_1_val <- seq(from=wr_sd_init[1]/2,to=wr_sd_init[1]*2,length.out=200)
wr_sd_lh_1_val <- wr_sd_1_val
wr_sd_vec <- wr_sd_init

for(i in 1:length(sigma_val)){
  wr_sd_vec[1] <- wr_sd_1_val[i]
  wr_sd_lh_1_val[i] <- log_likelihood_r(viral_data,
                                     individual_data,
                                     settings,
                                     wp_mean_init,
                                     wp_sd_init,
                                     dp_mean_init,
                                     dp_sd_init,
                                     tp_sd_init,
                                     wr_mean_init,
                                     wr_sd_vec,
                                     wp_init,
                                     dp_init,
                                     tp_init,
                                     wr_init,
                                     model_init,
                                     sigma_init)
}

plot(wr_sd_1_val,wr_sd_lh_1_val,type="l",main="WR sd 1 profile likelihood")
abline(v=wr_sd_init[1],col="blue")
print(wr_sd_1_val[which.max(wr_sd_lh_1_val)])

### WP Mean and SD
wp_mean_1_val <- seq(from=wp_mean_init[1]/2,to=wp_mean_init[1]*2,length.out=200)
wp_lh_1_val <- wp_mean_1_val
wp_mean_vec <- wp_mean_init

for(i in 1:length(sigma_val)){
  wp_mean_vec[1] <- wp_mean_1_val[i]
  wp_lh_1_val[i] <- log_likelihood_r(viral_data,
                                     individual_data,
                                     settings,
                                     wp_mean_vec,
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
                                     sigma_init)
}

plot(wp_mean_1_val,wp_lh_1_val,type="l",main="WP Mean 1 profile likelihood")
abline(v=wp_mean_init[1],col="blue")
print(wp_mean_1_val[which.max(wp_lh_1_val)])

wp_mean_2_val <- seq(from=wp_mean_init[2]/2,to=wp_mean_init[2]*2,length.out=200)
wp_lh_2_val <- wp_mean_2_val
wp_mean_vec <- wp_mean_init

for(i in 1:length(sigma_val)){
  wp_mean_vec[2] <- wp_mean_2_val[i]
  wp_lh_2_val[i] <- log_likelihood_r(viral_data,
                                     individual_data,
                                     settings,
                                     wp_mean_vec,
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
                                     sigma_init)
}

plot(wp_mean_2_val,wp_lh_2_val,type="l",main="WP Mean 2 profile likelihood")
abline(v=wp_mean_init[2],col="blue")
print(wp_mean_2_val[which.max(wp_lh_2_val)])

wp_mean_3_val <- seq(from=wp_mean_init[3]/2,to=wp_mean_init[3]*2,length.out=200)
wp_lh_3_val <- wp_mean_3_val
wp_mean_vec <- wp_mean_init

for(i in 1:length(sigma_val)){
  wp_mean_vec[3] <- wp_mean_3_val[i]
  wp_lh_3_val[i] <- log_likelihood_r(viral_data,
                                     individual_data,
                                     settings,
                                     wp_mean_vec,
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
                                     sigma_init)
}

plot(wp_mean_3_val,wp_lh_3_val,type="l",main="WP Mean 3 profile likelihood")
abline(v=wp_mean_init[3],col="blue")
print(wp_mean_3_val[which.max(wp_lh_3_val)])

wp_sd_1_val <- seq(from=wp_sd_init[1]/2,to=wp_sd_init[1]*2,length.out=200)
wp_sd_lh_1_val <- wp_sd_1_val
wp_sd_vec <- wp_sd_init

for(i in 1:length(sigma_val)){
  wp_sd_vec[1] <- wp_sd_1_val[i]
  wp_sd_lh_1_val[i] <- log_likelihood_r(viral_data,
                                        individual_data,
                                        settings,
                                        wp_mean_init,
                                        wp_sd_vec,
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
                                        sigma_init)
}

plot(wp_sd_1_val,wp_sd_lh_1_val,type="l",main="WP sd 1 profile likelihood")
abline(v=wp_sd_init[1],col="blue")
print(wp_sd_1_val[which.max(wp_sd_lh_1_val)])

### DP
dp_mean_1_val <- seq(from=dp_mean_init[1]/2,to=dp_mean_init[1]*2,length.out=200)
dp_lh_1_val <- dp_mean_1_val
dp_mean_vec <- dp_mean_init

for(i in 1:length(sigma_val)){
  dp_mean_vec[1] <- dp_mean_1_val[i]
  dp_lh_1_val[i] <- log_likelihood_r(viral_data,
                                     individual_data,
                                     settings,
                                     wp_mean_init,
                                     dp_sd_init,
                                     dp_mean_vec,
                                     dp_sd_init,
                                     tp_sd_init,
                                     wr_mean_init,
                                     wr_sd_init,
                                     wp_init,
                                     dp_init,
                                     tp_init,
                                     wr_init,
                                     model_init,
                                     sigma_init)
}

plot(dp_mean_1_val,dp_lh_1_val,type="l",main="DP Mean 1 profile likelihood")
abline(v=dp_mean_init[1],col="blue")
print(dp_mean_1_val[which.max(dp_lh_1_val)])

dp_mean_2_val <- seq(from=dp_mean_init[2]/2,to=dp_mean_init[2]*2,length.out=200)
dp_lh_2_val <- dp_mean_2_val
dp_mean_vec <- dp_mean_init

for(i in 1:length(sigma_val)){
  dp_mean_vec[2] <- dp_mean_2_val[i]
  dp_lh_2_val[i] <- log_likelihood_r(viral_data,
                                     individual_data,
                                     settings,
                                     wp_mean_init,
                                     dp_sd_init,
                                     dp_mean_vec,
                                     dp_sd_init,
                                     tp_sd_init,
                                     wr_mean_init,
                                     wr_sd_init,
                                     wp_init,
                                     dp_init,
                                     tp_init,
                                     wr_init,
                                     model_init,
                                     sigma_init)
}

plot(dp_mean_2_val,dp_lh_2_val,type="l",main="DP Mean 2 profile likelihood")
abline(v=dp_mean_init[2],col="blue")
print(dp_mean_2_val[which.max(dp_lh_2_val)])

dp_mean_3_val <- seq(from=dp_mean_init[3]/2,to=dp_mean_init[3]*2,length.out=200)
dp_lh_3_val <- dp_mean_3_val
dp_mean_vec <- dp_mean_init

for(i in 1:length(sigma_val)){
  dp_mean_vec[3] <- dp_mean_3_val[i]
  dp_lh_3_val[i] <- log_likelihood_r(viral_data,
                                     individual_data,
                                     settings,
                                     wp_mean_init,
                                     dp_sd_init,
                                     dp_mean_vec,
                                     dp_sd_init,
                                     tp_sd_init,
                                     wr_mean_init,
                                     wr_sd_init,
                                     wp_init,
                                     dp_init,
                                     tp_init,
                                     wr_init,
                                     model_init,
                                     sigma_init)
}

plot(dp_mean_3_val,dp_lh_3_val,type="l",main="DP Mean 3 profile likelihood")
abline(v=dp_mean_init[3],col="blue")
print(dp_mean_3_val[which.max(dp_lh_3_val)])

dp_sd_1_val <- seq(from=dp_sd_init[1]/2,to=dp_sd_init[1]*2,length.out=200)
dp_sd_lh_1_val <- dp_sd_1_val
dp_sd_vec <- dp_sd_init

for(i in 1:length(sigma_val)){
  dp_sd_vec[1] <- dp_sd_1_val[i]
  dp_sd_lh_1_val[i] <- log_likelihood_r(viral_data,
                                        individual_data,
                                        settings,
                                        wp_mean_init,
                                        wp_sd_init,
                                        dp_mean_init,
                                        dp_sd_vec,
                                        tp_sd_init,
                                        wr_mean_init,
                                        wr_sd_init,
                                        wp_init,
                                        dp_init,
                                        tp_init,
                                        wr_init,
                                        model_init,
                                        sigma_init)
}

plot(dp_sd_1_val,dp_sd_lh_1_val,type="l",main="DP sd 1 profile likelihood")
abline(v=dp_sd_init[1],col="blue")
print(dp_sd_1_val[which.max(dp_sd_lh_1_val)])

## TP
tp_sd_1_val <- seq(from=tp_sd_init[1]/2,to=tp_sd_init[1]*2,length.out=200)
tp_sd_lh_1_val <- tp_sd_1_val
tp_sd_vec <- tp_sd_init

for(i in 1:length(sigma_val)){
  tp_sd_vec[1] <- tp_sd_1_val[i]
  tp_sd_lh_1_val[i] <- log_likelihood_r(viral_data,
                                        individual_data,
                                        settings,
                                        wp_mean_init,
                                        wp_sd_init,
                                        dp_mean_init,
                                        dp_sd_init,
                                        tp_sd_vec,
                                        wr_mean_init,
                                        wr_sd_init,
                                        wp_init,
                                        dp_init,
                                        tp_init,
                                        wr_init,
                                        model_init,
                                        sigma_init)
}

plot(tp_sd_1_val,tp_sd_lh_1_val,type="l",main="TP sd 1 profile likelihood")
abline(v=tp_sd_init[1],col="blue")
print(tp_sd_1_val[which.max(tp_sd_lh_1_val)])


err <- sim_out$viral_load_data$viral_load_obs-sim_out$viral_load_data$viral_load_mu
err <- err[err!=0]

sqrt(sum(err^2)/(length(err)-1))

sum(sim_out$viral_load_data$viral_load_obs == 0 & sim_out$viral_load_data$viral_load_mu > 0)
err_trunc <- sim_out$viral_load_data$viral_load_mu[sim_out$viral_load_data$viral_load_obs == 0 & sim_out$viral_load_data$viral_load_mu > 0]

sigma_arg <- param_settings_in$sigma

sigma_args <- seq(from=0.5,to=10,length.out=200)
sigma_llh <- sigma_args

for(i in 1:length(sigma_args)){
  sigma_llh[i] <- sum(log(dnorm(err,mean=0,sd=sigma_args[i]))) # +
    # sum(log(1-pnorm(err_trunc,mean=0,sd=sigma_args[i])))
}

plot(sigma_args,sigma_llh,type="l",main="Sigma likelihood")
abline(v=sigma_arg,col="blue")
print(sigma_args[which.max(sigma_llh)])




