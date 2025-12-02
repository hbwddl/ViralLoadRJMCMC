#!/usr/bin/env Rscript
library(ViralLoadRJMCMC)
library(dplyr)
library(Rcpp)
library(truncnorm)
library(ggpubr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

csv_name <- args[1]
csv_line <- as.numeric(args[2])

# csv_name <- "~/Documents/Research/Within-Host/ViralLoadRJMCMC/test/batch_parameters.csv"
# csv_line <- 1

batch_params <- read.csv(csv_name,header=T)
run_params <- batch_params[csv_line,]

setwd("/projects/lau_projects/phylodynamics_hannah/Within_Host_Viral_Load/RJMCMC_Results")
# setwd("~/Documents/Research/Within-Host/RJMCMC_Results")

mcmc_seed <- run_params$seed
set.seed(mcmc_seed)

# analysis_dir <- paste0("~/Documents/Research/Within-Host/RJMCMC_Results/analysis_seed_",mcmc_seed)
analysis_dir <- paste0("/projects/lau_projects/phylodynamics_hannah/Within_Host_Viral_Load/RJMCMC_Results/analysis_seed_",mcmc_seed)
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

data_settings_in <- list(lod=run_params$lod,
                         n=run_params$n,
                         p_group=as.vector(as.numeric(unlist(strsplit(run_params$p_group,";")))),
                         t_obs=as.vector(as.numeric(unlist(strsplit(run_params$t_obs,";")))),
                         sensitivity=run_params$test_sensitivity)

param_settings_in <- list(p_model=as.vector(as.numeric(unlist(strsplit(run_params$p_model,";")))),
                          wp_mean=as.vector(as.numeric(unlist(strsplit(run_params$wp_mean,";")))),
                          wp_sd=as.vector(as.numeric(unlist(strsplit(run_params$wp_sd,";")))),
                          tp_sd=as.vector(as.numeric(unlist(strsplit(run_params$tp_sd,";")))),
                          dp_mean=as.vector(as.numeric(unlist(strsplit(run_params$dp_mean,";")))),
                          dp_sd=as.vector(as.numeric(unlist(strsplit(run_params$dp_sd,";")))),
                          wr_mean=as.vector(as.numeric(unlist(strsplit(run_params$wr_mean,";")))),
                          wr_sd=as.vector(as.numeric(unlist(strsplit(run_params$wr_sd,";")))),
                          sigma=run_params$sigma,
                          wp_min=run_params$wp_min,
                          wp_max=run_params$wp_max,
                          dp_min=run_params$dp_min,
                          tp_min=run_params$tp_min,
                          tp_max=run_params$tp_max,
                          wr_min=run_params$wr_min,
                          wr_max=run_params$wr_max)

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

settings <- data.frame(lod=run_params$lod,
                       sensitivity=run_params$test_sensitivity,
                       n_iterations=run_params$n_iterations,
                       n_subtypes=n_subtype,
                       n_subjects=nrow(individual_data),
                       n_data=nrow(viral_data))

priors <- data.frame(wp_min = run_params$wp_min_prior,
                     wp_max = run_params$wp_max_prior,
                     wr_min = run_params$wr_min_prior,
                     wr_max = run_params$wr_max_prior,
                     wpmean_max = run_params$wpmean_max,
                     dpmean_max = run_params$dpmean_max,
                     wrmean_max = run_params$wrmean_max,
                     wpsd_max = run_params$wpsd_max,
                     tpsd_max = run_params$tpsd_max,
                     dpsd_max = run_params$dpsd_max,
                     wrsd_max = run_params$wrsd_max,
                     sigma_max = run_params$sigma_max,
                     wpsd_min = run_params$wpsd_min,
                     tpsd_min = run_params$tpsd_min,
                     dpsd_min = run_params$dpsd_min,
                     wrsd_min = run_params$wrsd_min,
                     sigma_min = run_params$sigma_min,
                     wpmean_mean = run_params$wpmean_mean,
                     wpmean_sd = run_params$wpmean_sd,
                     dpmean_mean = run_params$dpmean_mean,
                     dpmean_sd = run_params$dpmean_sd,
                     wrmean_mean = run_params$wrmean_mean,
                     wrmean_sd = run_params$wrmean_sd,
                     wpsd_scale = run_params$wpsd_scale,
                     tpsd_scale = run_params$tpsd_scale,
                     dpsd_scale = run_params$dpsd_scale,
                     wrsd_scale = run_params$wrsd_scale,
                     sigma_scale = run_params$sigma_scale
                     )

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
model_init <- rep(2,nrow(individual_data))
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

time1 <- Sys.time()

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

print(paste0("Total time: ",Sys.time() - time1))

sink(file=NULL)

quantile_coverage <- function(vec,true_val){
  return(quantile(vec,probs=c(0.025)) <= true_val & quantile(vec,probs=c(0.975)) >= true_val)
}

dec.precision <- 6

pct_burnin_begin <- 0.4
pct_burnin_end <- 0.99

load("individual_data_in.RData")
load("priors.RData")
load("settings_in.RData")
load("viral_data_in.RData")
load("param_settings_in.RData")

index_id <- 0:(nrow(individual_data)-1)

scalars_out_raw <- read.csv("./output/scalars_out.csv",header=T)
n_burnin_begin <- round(nrow(scalars_out_raw)*pct_burnin_begin)
n_burnin_end <- round(nrow(scalars_out_raw)*pct_burnin_end)

scalars_out <- scalars_out_raw[n_burnin_begin:n_burnin_end,]
rm(scalars_out_raw)

scalar_plotnames <- c("wp Mean, H1N1",
                      "wp Mean, H3N2",
                      "wp Mean, Dual",
                      "wp SD, H1N1",
                      "wp SD, H3N2",
                      "wp SD, Dual",
                      "dp Mean, H1N1",
                      "dp Mean, H3N2",
                      "dp Mean, Dual",
                      "dp SD, H1N1",
                      "dp SD, H3N2",
                      "dp SD, Dual",
                      "tp SD, H1N1",
                      "tp SD, H3N2",
                      "tp SD, Dual",
                      "wr Mean, H1N1",
                      "wr Mean, H3N2",
                      "wr Mean, Dual",
                      "wr SD, H1N1",
                      "wr SD, H3N2",
                      "wr SD, Dual",
                      "Sigma",
                      "Log Likelihood")

print("Scalar plots")

pdf(file="Scalar_Plots.pdf",width=12,height=8)
par(mfrow=c(2,3))

for(i in 1:ncol(scalars_out)){
  print(scalar_plotnames[i])
  plot(scalars_out[,i],type="l",main=scalar_plotnames[i])
}

par(mfrow=c(1,1))
dev.off()

print("Scalar Difference Plots")

pdf(file="Scalar_Diff_Plots.pdf",width=14,height=5)

par(mfrow=c(1,3))
plot(scalars_out$wp_mean_0 - scalars_out$wp_mean_1,type="l",main="WP Mean Difference, H1N1-H3N2",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$wp_mean_2 - scalars_out$wp_mean_0,type="l",main="WP Mean Difference, Dual-H1N1",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$wp_mean_2 - scalars_out$wp_mean_1,type="l",main="WP Mean Difference, Dual-H3N2",xlab="Iteration",ylab="Difference")
abline(h=0)

plot(scalars_out$dp_mean_0 - scalars_out$dp_mean_1,type="l",main="DP Mean Difference, H1N1-H3N2",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$dp_mean_2 - scalars_out$dp_mean_0,type="l",main="DP Mean Difference, Dual-H1N1",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$dp_mean_2 - scalars_out$dp_mean_1,type="l",main="DP Mean Difference, Dual-H3N2",xlab="Iteration",ylab="Difference")
abline(h=0)

plot(scalars_out$wr_mean_0 - scalars_out$wr_mean_1,type="l",main="WR Mean Difference, H1N1-H3N2",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$wr_mean_2 - scalars_out$wr_mean_0,type="l",main="WR Mean Difference, Dual-H1N1",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$wr_mean_2 - scalars_out$wr_mean_1,type="l",main="WR Mean Difference, Dual-H3N2",xlab="Iteration",ylab="Difference")
abline(h=0)
par(mfrow=c(1,1))

dev.off()

print("Quantiles")

sink(file="Quantiles.txt")
print("wp Mean, h1n1")
print(quantile(scalars_out$wp_mean_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_mean_0,param_settings_in$wp_mean[1]))

print("wp Mean, h3n2")
print(quantile(scalars_out$wp_mean_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_mean_1,param_settings_in$wp_mean[2]))

print("wp Mean, dual")
print(quantile(scalars_out$wp_mean_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_mean_2,param_settings_in$wp_mean[3]))

print("dp Mean, h1n1")
print(quantile(scalars_out$dp_mean_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_mean_0,param_settings_in$dp_mean[1]))

print("dp Mean, h3n2")
print(quantile(scalars_out$dp_mean_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_mean_1,param_settings_in$dp_mean[2]))

print("dp Mean, dual")
print(quantile(scalars_out$dp_mean_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_mean_2,param_settings_in$dp_mean[3]))

print("wr Mean, h1n1")
print(quantile(scalars_out$wr_mean_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_mean_0,param_settings_in$wr_mean[1]))

print("wr Mean, h3n2")
print(quantile(scalars_out$wr_mean_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_mean_1,param_settings_in$wr_mean[2]))

print("wr Mean, dual")
print(quantile(scalars_out$wr_mean_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_mean_2,param_settings_in$wr_mean[3]))

print("WP Mean Difference, H1N1-H3N2")
print(quantile(scalars_out$wp_mean_0 - scalars_out$wp_mean_1,probs=c(0.025,0.5,0.975)))

print("WP Mean Difference, Dual-H1N1")
print(quantile(scalars_out$wp_mean_2 - scalars_out$wp_mean_0,probs=c(0.025,0.5,0.975)))

print("WP Mean Difference, Dual-H3N2")
print(quantile(scalars_out$wp_mean_2 - scalars_out$wp_mean_1,probs=c(0.025,0.5,0.975)))

print("DP Mean Difference, H1N1-H3N2")
print(quantile(scalars_out$dp_mean_0 - scalars_out$dp_mean_1,probs=c(0.025,0.5,0.975)))

print("DP Mean Difference, Dual-H1N1")
print(quantile(scalars_out$dp_mean_2 - scalars_out$dp_mean_0,probs=c(0.025,0.5,0.975)))

print("DP Mean Difference, Dual-H3N2")
print(quantile(scalars_out$dp_mean_2 - scalars_out$dp_mean_1,probs=c(0.025,0.5,0.975)))

print("WR Mean Difference, H1N1-H3N2")
print(quantile(scalars_out$wr_mean_0 - scalars_out$wr_mean_1,probs=c(0.025,0.5,0.975)))

print("WR Mean Difference, Dual-H1N1")
print(quantile(scalars_out$wr_mean_2 - scalars_out$wr_mean_0,probs=c(0.025,0.5,0.975)))

print("WR Mean Difference, Dual-H3N2")
print(quantile(scalars_out$wr_mean_2 - scalars_out$wr_mean_1,probs=c(0.025,0.5,0.975)))

print("wp SD, h1n1")
print(quantile(scalars_out$wp_sd_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_sd_0,param_settings_in$wp_sd[1]))

print("wp SD, h3n2")
print(quantile(scalars_out$wp_sd_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_sd_1,param_settings_in$wp_sd[2]))

print("wp SD, dual")
print(quantile(scalars_out$wp_sd_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_sd_2,param_settings_in$wp_sd[3]))

print("tp SD, h1n1")
print(quantile(scalars_out$tp_sd_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$tp_sd_0,param_settings_in$tp_sd[1]))

print("tp SD, h3n2")
print(quantile(scalars_out$tp_sd_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$tp_sd_1,param_settings_in$tp_sd[2]))

print("tp SD, dual")
print(quantile(scalars_out$tp_sd_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$tp_sd_2,param_settings_in$tp_sd[3]))

print("dp SD, h1n1")
print(quantile(scalars_out$dp_sd_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_sd_0,param_settings_in$dp_sd[1]))

print("dp SD, h3n2")
print(quantile(scalars_out$dp_sd_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_sd_1,param_settings_in$dp_sd[2]))

print("dp SD, dual")
print(quantile(scalars_out$dp_sd_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_sd_2,param_settings_in$dp_sd[3]))

print("wr SD, h1n1")
print(quantile(scalars_out$wr_sd_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_sd_0,param_settings_in$wr_sd[1]))

print("wr SD, h3n2")
print(quantile(scalars_out$wr_sd_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_sd_1,param_settings_in$wr_sd[2]))

print("wr SD, dual")
print(quantile(scalars_out$wr_sd_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_sd_2,param_settings_in$wr_sd[3]))

print("sigma")
print(quantile(scalars_out$sigma,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$sigma,param_settings_in$sigma))

acp_pr <- function(mcmc_vec){
  return(mean(mcmc_vec[1:(length(mcmc_vec)-1)] != mcmc_vec[2:(length(mcmc_vec))],na.rm=T))
}

scalar_acp_pr <- apply(scalars_out,2,acp_pr)

print("Acceptance Probabilities")
print(scalar_acp_pr)

sink(file=NULL)

model_out_raw <- read.csv("./output/model_out.csv",header=F)
n_burnin_begin <- round(nrow(model_out_raw)*pct_burnin_begin)
n_burnin_end <- round(nrow(model_out_raw)*pct_burnin_end)

model_out <- model_out_raw[n_burnin_begin:n_burnin_end,]
rm(model_out_raw)

wp_out_raw <- read.csv("./output/wp_out.csv",header=F)
n_burnin_begin <- round(nrow(wp_out_raw)*pct_burnin_begin)
n_burnin_end <- round(nrow(wp_out_raw)*pct_burnin_end)

wp_out <- wp_out_raw[n_burnin_begin:n_burnin_end,]
rm(wp_out_raw)

tp_out_raw <- read.csv("./output/tp_out.csv",header=F)
n_burnin_begin <- round(nrow(tp_out_raw)*pct_burnin_begin)
n_burnin_end <- round(nrow(tp_out_raw)*pct_burnin_end)

tp_out <- tp_out_raw[n_burnin_begin:n_burnin_end,]
rm(tp_out_raw)

dp_out_raw <- read.csv("./output/dp_out.csv",header=F)
n_burnin_begin <- round(nrow(dp_out_raw)*pct_burnin_begin)
n_burnin_end <- round(nrow(dp_out_raw)*pct_burnin_end)

dp_out <- dp_out_raw[n_burnin_begin:n_burnin_end,]
rm(dp_out_raw)

wr_out_raw <- read.csv("./output/wr_out.csv",header=F)
n_burnin_begin <- round(nrow(wr_out_raw)*pct_burnin_begin)
n_burnin_end <- round(nrow(wr_out_raw)*pct_burnin_end)

wr_out <- wr_out_raw[n_burnin_begin:n_burnin_end,]
rm(wr_out_raw)

print("Model Traceplots")

pdf(file="Model_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))

for(i in 1:ncol(model_out)){
  model_p <- round(table(c(model_out[,i],1,2,3))/(length(model_out[,i])+3),3)
  
  plot(model_out[,i],main=paste0("Model, ",index_id[i],", p=(",model_p[1],",",model_p[2],",",model_p[3],")"),type="l",xlab="Iteration",ylab="Model")
}

par(mfrow=c(1,1))

dev.off()

print("WP Traceplots")

pdf(file="WP_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))
for(i in 1:ncol(wp_out)){
  plot(wp_out[,i],main=paste0("wp, individual ",index_id[i]),type="l",xlab="Iteration",ylab="wp")
}
par(mfrow=c(1,1))

dev.off()

print("TP Traceplots")

pdf(file="TP_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))
for(i in 1:ncol(tp_out)){
  plot(tp_out[,i],main=paste0("Tp, individual ",index_id[i]),type="l",xlab="Iteration",ylab="tp")
}
par(mfrow=c(1,1))

dev.off()

print("DP Traceplots")

pdf(file="DP_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))
for(i in 1:ncol(dp_out)){
  plot(dp_out[,i],main=paste0("dp, individual ",index_id[i]),type="l",xlab="Iteration",ylab="dp")
}
par(mfrow=c(1,1))

dev.off()

print("WR Traceplots")

pdf(file="WR_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))
for(i in 1:ncol(wr_out)){
  plot(wr_out[,i],main=paste0("wr, individual ",index_id[i]),type="l",xlab="Iteration",ylab="wr")
}
par(mfrow=c(1,1))

dev.off()

print("Data Plots")

pdf(file="Data_plots.pdf",width=16,height=16)
par(mfrow=c(5,5))

for(i in 1:nrow(individual_data)){
  plot_dat <- viral_data %>%
    filter(index_init==individual_data$index_init[i]) %>%
    arrange(time)
  
  wp_quantile <- quantile(wp_out[,i],probs = c(0.025,0.5,0.975))
  tp_quantile <- quantile(tp_out[,i],probs = c(0.025,0.5,0.975))
  dp_quantile <- quantile(dp_out[,i],probs = c(0.025,0.5,0.975))
  wr_quantile <- quantile(wr_out[,i],probs = c(0.025,0.5,0.975))
  
  tp_med <- median(tp_out[,i])
  
  model_infer <- which.max(table(c(model_out[,i],1,2,3)))
  
  plot_col <- c("blue","red","purple")[individual_data$subtype[i]+1]
  
  # plot_col <- "blue"
  
  plot(plot_dat$time,plot_dat$viral_load,pch=19,main=paste0("Observed Data, ID ",i),col=plot_col,
       xlim=c(min(-wp_quantile[3],plot_dat$time),max(wr_quantile[3],plot_dat$time)),
       ylim=c(0,max(dp_quantile[3],plot_dat$time)))
  
  if(model_infer != 3){
    lines(c(-wp_quantile[2],tp_med),c(0,dp_quantile[2]))
    lines(c(-wp_quantile[1],tp_med),c(0,dp_quantile[1]),lty="dashed")
    lines(c(-wp_quantile[3],tp_med),c(0,dp_quantile[3]),lty="dashed")
  }
  if(model_infer != 1){
    lines(c(tp_med,wr_quantile[2]),c(dp_quantile[2],0))
    lines(c(tp_med,wr_quantile[1]),c(dp_quantile[1],0),lty="dashed")
    lines(c(tp_med,wr_quantile[3]),c(dp_quantile[3],0),lty="dashed")
  }
  
  
}

par(mfrow=c(1,1))
dev.off()

quantile_coverage <- data.frame(wp_mean_0_lo=quantile(scalars_out$wp_mean_0,probs=c(0.025)),
                                wp_mean_0_hi=quantile(scalars_out$wp_mean_0,probs=c(0.975)),
                                wp_mean_0_cover=quantile_coverage(scalars_out$wp_mean_0,param_settings_in$wp_mean[1]),
                                wp_mean_1_lo=quantile(scalars_out$wp_mean_1,probs=c(0.025)),
                                wp_mean_1_hi=quantile(scalars_out$wp_mean_1,probs=c(0.975)),
                                wp_mean_1_cover=quantile_coverage(scalars_out$wp_mean_1,param_settings_in$wp_mean[2]),
                                wp_mean_2_lo=quantile(scalars_out$wp_mean_2,probs=c(0.025)),
                                wp_mean_2_hi=quantile(scalars_out$wp_mean_2,probs=c(0.975)),
                                wp_mean_2_cover=quantile_coverage(scalars_out$wp_mean_2,param_settings_in$wp_mean[3]),
                                dp_mean_0_lo=quantile(scalars_out$dp_mean_0,probs=c(0.025)),
                                dp_mean_0_hi=quantile(scalars_out$dp_mean_0,probs=c(0.975)),
                                dp_mean_0_cover=quantile_coverage(scalars_out$dp_mean_0,param_settings_in$dp_mean[1]),
                                dp_mean_1_lo=quantile(scalars_out$dp_mean_1,probs=c(0.025)),
                                dp_mean_1_hi=quantile(scalars_out$dp_mean_1,probs=c(0.975)),
                                dp_mean_1_cover=quantile_coverage(scalars_out$dp_mean_1,param_settings_in$dp_mean[2]),
                                dp_mean_2_lo=quantile(scalars_out$dp_mean_2,probs=c(0.025)),
                                dp_mean_2_hi=quantile(scalars_out$dp_mean_2,probs=c(0.975)),
                                dp_mean_2_cover=quantile_coverage(scalars_out$dp_mean_2,param_settings_in$dp_mean[3]),
                                wr_mean_0_lo=quantile(scalars_out$wr_mean_0,probs=c(0.025)),
                                wr_mean_0_hi=quantile(scalars_out$wr_mean_0,probs=c(0.975)),
                                wr_mean_0_cover=quantile_coverage(scalars_out$wr_mean_0,param_settings_in$wr_mean[1]),
                                wr_mean_1_lo=quantile(scalars_out$wr_mean_1,probs=c(0.025)),
                                wr_mean_1_hi=quantile(scalars_out$wr_mean_1,probs=c(0.975)),
                                wr_mean_1_cover=quantile_coverage(scalars_out$wr_mean_1,param_settings_in$wr_mean[2]),
                                wr_mean_2_lo=quantile(scalars_out$wr_mean_2,probs=c(0.025)),
                                wr_mean_2_hi=quantile(scalars_out$wr_mean_2,probs=c(0.975)),
                                wr_mean_2_cover=quantile_coverage(scalars_out$wr_mean_2,param_settings_in$wr_mean[3]),
                                wp_sd_0_lo=quantile(scalars_out$wp_sd_0,probs=c(0.025)),
                                wp_sd_0_hi=quantile(scalars_out$wp_sd_0,probs=c(0.975)),
                                wp_sd_0_cover=quantile_coverage(scalars_out$wp_sd_0,param_settings_in$wp_sd[1]),
                                wp_sd_1_lo=quantile(scalars_out$wp_sd_1,probs=c(0.025)),
                                wp_sd_1_hi=quantile(scalars_out$wp_sd_1,probs=c(0.975)),
                                wp_sd_1_cover=quantile_coverage(scalars_out$wp_sd_1,param_settings_in$wp_sd[2]),
                                wp_sd_2_lo=quantile(scalars_out$wp_sd_2,probs=c(0.025)),
                                wp_sd_2_hi=quantile(scalars_out$wp_sd_2,probs=c(0.975)),
                                wp_sd_2_cover=quantile_coverage(scalars_out$wp_sd_2,param_settings_in$wp_mean[3]),
                                tp_sd_0_lo=quantile(scalars_out$tp_sd_0,probs=c(0.025)),
                                tp_sd_0_hi=quantile(scalars_out$tp_sd_0,probs=c(0.975)),
                                tp_sd_0_cover=quantile_coverage(scalars_out$tp_sd_0,param_settings_in$tp_sd[1]),
                                tp_sd_1_lo=quantile(scalars_out$tp_sd_1,probs=c(0.025)),
                                tp_sd_1_hi=quantile(scalars_out$tp_sd_1,probs=c(0.975)),
                                tp_sd_1_cover=quantile_coverage(scalars_out$tp_sd_1,param_settings_in$tp_sd[2]),
                                tp_sd_2_lo=quantile(scalars_out$tp_sd_2,probs=c(0.025)),
                                tp_sd_2_hi=quantile(scalars_out$tp_sd_2,probs=c(0.975)),
                                tp_sd_2_cover=quantile_coverage(scalars_out$tp_sd_2,param_settings_in$tp_sd[3]),
                                dp_sd_0_lo=quantile(scalars_out$dp_sd_0,probs=c(0.025)),
                                dp_sd_0_hi=quantile(scalars_out$dp_sd_0,probs=c(0.975)),
                                dp_sd_0_cover=quantile_coverage(scalars_out$dp_sd_0,param_settings_in$dp_sd[1]),
                                dp_sd_1_lo=quantile(scalars_out$dp_sd_1,probs=c(0.025)),
                                dp_sd_1_hi=quantile(scalars_out$dp_sd_1,probs=c(0.975)),
                                dp_sd_1_cover=quantile_coverage(scalars_out$dp_sd_1,param_settings_in$dp_sd[2]),
                                dp_sd_2_lo=quantile(scalars_out$dp_sd_2,probs=c(0.025)),
                                dp_sd_2_hi=quantile(scalars_out$dp_sd_2,probs=c(0.975)),
                                dp_sd_2_cover=quantile_coverage(scalars_out$dp_sd_2,param_settings_in$dp_sd[3]),
                                wr_sd_0_lo=quantile(scalars_out$wr_sd_0,probs=c(0.025)),
                                wr_sd_0_hi=quantile(scalars_out$wr_sd_0,probs=c(0.975)),
                                wr_sd_0_cover=quantile_coverage(scalars_out$wr_sd_0,param_settings_in$wr_sd[1]),
                                wr_sd_1_lo=quantile(scalars_out$wr_sd_1,probs=c(0.025)),
                                wr_sd_1_hi=quantile(scalars_out$wr_sd_1,probs=c(0.975)),
                                wr_sd_1_cover=quantile_coverage(scalars_out$wr_sd_1,param_settings_in$wr_sd[2]),
                                wr_sd_2_lo=quantile(scalars_out$wr_sd_2,probs=c(0.025)),
                                wr_sd_2_hi=quantile(scalars_out$wr_sd_2,probs=c(0.975)),
                                wr_sd_2_cover=quantile_coverage(scalars_out$wr_sd_2,param_settings_in$wr_sd[3]),
                                sigma_lo=quantile(scalars_out$sigma,probs=c(0.025)),
                                sigma_hi=quantile(scalars_out$sigma,probs=c(0.975)),
                                sigma_cover=quantile_coverage(scalars_out$sigma,param_settings_in$sigma))

save(quantile_coverage,file="quantile_coverage.RData")