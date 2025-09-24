library(Rcpp)
library(dplyr)
library(truncnorm)
library(ggpubr)
library(ggplot2)

quantile_coverage <- function(vec,true_val){
  return(quantile(vec,probs=c(0.025)) <= true_val & quantile(vec,probs=c(0.975)) >= true_val)
}

dec.precision <- 6

pct_burnin_begin <- 0.0
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
  
  plot_col <- c("blue","red","purple")[individual_data$subtype+1]
  
  plot_col <- "blue"
  
  plot(plot_dat$time,plot_dat$viral_load,pch=19,main=paste0("Observed Data, ID ",individual_data$index_init[i]),col=plot_col,
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

# print("Model counts")
# 
# model_counts <- matrix(data=NA,nrow=nrow(model_out),ncol=3)
# 
# for(i in 1:nrow(model_out)){
#   model_p <- table(unlist(model_out[i,]))
#   
#   model_counts[i,1] <- model_p[1]
#   model_counts[i,2] <- model_p[2]
#   model_counts[i,3] <- model_p[3]
# }
# 
# model_counts[is.na(model_counts)] <- 0
# 
# plot(model_counts[,1],type="l",ylim=c(0,max(model_counts)),main="# of individuals in each model")
# lines(model_counts[,2],col="blue")
# lines(model_counts[,3],col="red")
# 
