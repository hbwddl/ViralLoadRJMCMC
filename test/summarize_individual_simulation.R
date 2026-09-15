library(Rcpp)
library(dplyr)
library(truncnorm)
library(ggpubr)
library(ggplot2)
library(coda)
library(ViralLoadRJMCMC)

setwd("~/Documents/Research/Within-Host/RJMCMC_Results/model_switch_results_good/analysis_seed_1")
# setwd("~/Documents/Research/Within-Host/RJMCMC_Results/analysis_seed_1")

quantile_coverage <- function(vec,true_val){
  return(quantile(vec,probs=c(0.025)) <= true_val & quantile(vec,probs=c(0.975)) >= true_val)
}

dec.precision <- 6

seeds <- c(1)

pct_burnin_begin <- 0.5
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

scalar_plotnames <- c("wp Mean, 0",
                      "wp Mean, 1",
                      "wp Mean, 2",
                      "wp SD, 0",
                      "wp SD, 1",
                      "wp SD, 2",
                      "dp Mean, 0",
                      "dp Mean, 1",
                      "dp Mean, 2",
                      "dp SD, 0",
                      "dp SD, 1",
                      "dp SD, 2",
                      "tp SD, 0",
                      "tp SD, 1",
                      "tp SD, 2",
                      "wr Mean, 0",
                      "wr Mean, 1",
                      "wr Mean, 2",
                      "wr SD, 0",
                      "wr SD, 1",
                      "wr SD, 2",
                      "Sigma",
                      "Log Likelihood")

scalar_true_values <- c(param_settings_in$wp_mean[1],
                        param_settings_in$wp_mean[2],
                        param_settings_in$wp_mean[3],
                        param_settings_in$wp_sd[1],
                        param_settings_in$wp_sd[2],
                        param_settings_in$wp_sd[3],
                        param_settings_in$dp_mean[1],
                        param_settings_in$dp_mean[2],
                        param_settings_in$dp_mean[3],
                        param_settings_in$dp_sd[1],
                        param_settings_in$dp_sd[2],
                        param_settings_in$dp_sd[3],
                        param_settings_in$tp_sd[1],
                        param_settings_in$tp_sd[2],
                        param_settings_in$tp_sd[3],
                        param_settings_in$wr_mean[1],
                        param_settings_in$wr_mean[2],
                        param_settings_in$wr_mean[3],
                        param_settings_in$wr_sd[1],
                        param_settings_in$wr_sd[2],
                        param_settings_in$wr_sd[3],
                        param_settings_in$sigma,
                        NA)

print("Scalar plots")

pdf(file="Scalar_Plots.pdf",width=12,height=8)
par(mfrow=c(2,3))

for(i in 2:ncol(scalars_out)){
  print(scalar_plotnames[i-1])
  plot(scalars_out[,i],type="l",main=scalar_plotnames[i-1])
  abline(h=scalar_true_values[i-1],col="blue")
}

par(mfrow=c(1,1))
dev.off()

print("Scalar Difference Plots")

pdf(file="Scalar_Diff_Plots.pdf",width=14,height=5)

par(mfrow=c(1,3))
plot(scalars_out$wp_mean_0 - scalars_out$wp_mean_1,type="l",main="WP Mean Difference, Subtype 1-Subtype 2",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$wp_mean_2 - scalars_out$wp_mean_0,type="l",main="WP Mean Difference, Subtype 3-Subtype 1",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$wp_mean_2 - scalars_out$wp_mean_1,type="l",main="WP Mean Difference, Subtype 3-Subtype 2",xlab="Iteration",ylab="Difference")
abline(h=0)

plot(scalars_out$dp_mean_0 - scalars_out$dp_mean_1,type="l",main="DP Mean Difference, Subtype 1-Subtype 2",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$dp_mean_2 - scalars_out$dp_mean_0,type="l",main="DP Mean Difference, Subtype 3-Subtype 1",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$dp_mean_2 - scalars_out$dp_mean_1,type="l",main="DP Mean Difference, Subtype 3-Subtype 2",xlab="Iteration",ylab="Difference")
abline(h=0)

plot(scalars_out$wr_mean_0 - scalars_out$wr_mean_1,type="l",main="WR Mean Difference, Subtype 1-Subtype 2",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$wr_mean_2 - scalars_out$wr_mean_0,type="l",main="WR Mean Difference, Subtype 3-Subtype 1",xlab="Iteration",ylab="Difference")
abline(h=0)
plot(scalars_out$wr_mean_2 - scalars_out$wr_mean_1,type="l",main="WR Mean Difference, Subtype 3-Subtype 2",xlab="Iteration",ylab="Difference")
abline(h=0)
par(mfrow=c(1,1))

dev.off()

quantile_cover <- data.frame("wp_mean_0" = NA,
                             "wp_mean_1" = NA,
                             "wp_mean_2" = NA,
                             "wp_sd_0" = NA,
                             "wp_sd_1" = NA,
                             "wp_sd_2" = NA,
                             "dp_mean_0" = NA,
                             "dp_mean_1" = NA,
                             "dp_mean_2" = NA,
                             "dp_sd_0" = NA,
                             "dp_sd_1" = NA,
                             "dp_sd_2" = NA,
                             "wr_mean_0" = NA,
                             "wr_mean_1" = NA,
                             "wr_mean_2" = NA,
                             "wr_sd_0" = NA,
                             "wr_sd_1" = NA,
                             "wr_sd_2" = NA,
                             "sigma" = NA)

print("Quantiles")

sink(file="Quantiles.txt")
print("wp Mean, h1n1")
print(quantile(scalars_out$wp_mean_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_mean_0,param_settings_in$wp_mean[1]))
quantile_cover$wp_mean_0 = quantile_coverage(scalars_out$wp_mean_0,param_settings_in$wp_mean[1])

print("wp Mean, h3n2")
print(quantile(scalars_out$wp_mean_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_mean_1,param_settings_in$wp_mean[2]))
quantile_cover$wp_mean_1 = quantile_coverage(scalars_out$wp_mean_1,param_settings_in$wp_mean[2])

print("wp Mean, dual")
print(quantile(scalars_out$wp_mean_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_mean_2,param_settings_in$wp_mean[3]))
quantile_cover$wp_mean_2 = quantile_coverage(scalars_out$wp_mean_2,param_settings_in$wp_mean[3])

print("dp Mean, h1n1")
print(quantile(scalars_out$dp_mean_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_mean_0,param_settings_in$dp_mean[1]))
quantile_cover$dp_mean_0 = quantile_coverage(scalars_out$dp_mean_0,param_settings_in$dp_mean[1])

print("dp Mean, h3n2")
print(quantile(scalars_out$dp_mean_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_mean_1,param_settings_in$dp_mean[2]))
quantile_cover$dp_mean_1 = quantile_coverage(scalars_out$dp_mean_1,param_settings_in$dp_mean[2])

print("dp Mean, dual")
print(quantile(scalars_out$dp_mean_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_mean_2,param_settings_in$dp_mean[3]))
quantile_cover$dp_mean_2 = quantile_coverage(scalars_out$dp_mean_2,param_settings_in$dp_mean[3])

print("wr Mean, h1n1")
print(quantile(scalars_out$wr_mean_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_mean_0,param_settings_in$wr_mean[1]))
quantile_cover$wr_mean_0 = quantile_coverage(scalars_out$wr_mean_0,param_settings_in$wr_mean[1])

print("wr Mean, h3n2")
print(quantile(scalars_out$wr_mean_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_mean_1,param_settings_in$wr_mean[2]))
quantile_cover$wr_mean_1 = quantile_coverage(scalars_out$wr_mean_1,param_settings_in$wr_mean[2])

print("wr Mean, dual")
print(quantile(scalars_out$wr_mean_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_mean_2,param_settings_in$wr_mean[3]))
quantile_cover$wr_mean_2 = quantile_coverage(scalars_out$wr_mean_2,param_settings_in$wr_mean[3])

print("WP Mean Difference, Subtype 1-Subtype 2")
print(quantile(scalars_out$wp_mean_0 - scalars_out$wp_mean_1,probs=c(0.025,0.5,0.975)))

print("WP Mean Difference, Subtype 3-Subtype 1")
print(quantile(scalars_out$wp_mean_2 - scalars_out$wp_mean_0,probs=c(0.025,0.5,0.975)))

print("WP Mean Difference, Subtype 3-Subtype 2")
print(quantile(scalars_out$wp_mean_2 - scalars_out$wp_mean_1,probs=c(0.025,0.5,0.975)))

print("DP Mean Difference, Subtype 1-Subtype 2")
print(quantile(scalars_out$dp_mean_0 - scalars_out$dp_mean_1,probs=c(0.025,0.5,0.975)))

print("DP Mean Difference, Subtype 3-Subtype 1")
print(quantile(scalars_out$dp_mean_2 - scalars_out$dp_mean_0,probs=c(0.025,0.5,0.975)))

print("DP Mean Difference, Subtype 3-Subtype 2")
print(quantile(scalars_out$dp_mean_2 - scalars_out$dp_mean_1,probs=c(0.025,0.5,0.975)))

print("WR Mean Difference, Subtype 1-Subtype 2")
print(quantile(scalars_out$wr_mean_0 - scalars_out$wr_mean_1,probs=c(0.025,0.5,0.975)))

print("WR Mean Difference, Subtype 3-Subtype 1")
print(quantile(scalars_out$wr_mean_2 - scalars_out$wr_mean_0,probs=c(0.025,0.5,0.975)))

print("WR Mean Difference, Subtype 3-Subtype 2")
print(quantile(scalars_out$wr_mean_2 - scalars_out$wr_mean_1,probs=c(0.025,0.5,0.975)))

print("wp SD, h1n1")
print(quantile(scalars_out$wp_sd_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_sd_0,param_settings_in$wp_sd[1]))
quantile_cover$wp_sd_0 = quantile_coverage(scalars_out$wp_sd_0,param_settings_in$wp_sd[1])

print("wp SD, h3n2")
print(quantile(scalars_out$wp_sd_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_sd_1,param_settings_in$wp_sd[2]))
quantile_cover$wp_sd_1 = quantile_coverage(scalars_out$wp_sd_1,param_settings_in$wp_sd[2])

print("wp SD, dual")
print(quantile(scalars_out$wp_sd_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wp_sd_2,param_settings_in$wp_sd[3]))
quantile_cover$wp_sd_2 = quantile_coverage(scalars_out$wp_sd_2,param_settings_in$wp_sd[3])

print("tp SD, h1n1")
print(quantile(scalars_out$tp_sd_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$tp_sd_0,param_settings_in$tp_sd[1]))
quantile_cover$tp_sd_0 = quantile_coverage(scalars_out$tp_sd_0,param_settings_in$tp_sd[1])

print("tp SD, h3n2")
print(quantile(scalars_out$tp_sd_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$tp_sd_1,param_settings_in$tp_sd[2]))
quantile_cover$tp_sd_1 = quantile_coverage(scalars_out$tp_sd_1,param_settings_in$tp_sd[2])

print("tp SD, dual")
print(quantile(scalars_out$tp_sd_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$tp_sd_2,param_settings_in$tp_sd[3]))
quantile_cover$tp_sd_2 = quantile_coverage(scalars_out$tp_sd_2,param_settings_in$tp_sd[3])

print("dp SD, h1n1")
print(quantile(scalars_out$dp_sd_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_sd_0,param_settings_in$dp_sd[1]))
quantile_cover$dp_sd_0 = quantile_coverage(scalars_out$dp_sd_0,param_settings_in$dp_sd[1])

print("dp SD, h3n2")
print(quantile(scalars_out$dp_sd_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_sd_1,param_settings_in$dp_sd[2]))
quantile_cover$dp_sd_1 = quantile_coverage(scalars_out$dp_sd_1,param_settings_in$dp_sd[2])

print("dp SD, dual")
print(quantile(scalars_out$dp_sd_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$dp_sd_2,param_settings_in$dp_sd[3]))
quantile_cover$dp_sd_2 = quantile_coverage(scalars_out$dp_sd_2,param_settings_in$dp_sd[3])

print("wr SD, h1n1")
print(quantile(scalars_out$wr_sd_0,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_sd_0,param_settings_in$wr_sd[1]))
quantile_cover$wr_sd_0 = quantile_coverage(scalars_out$wr_sd_0,param_settings_in$wr_sd[1])

print("wr SD, h3n2")
print(quantile(scalars_out$wr_sd_1,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_sd_1,param_settings_in$wr_sd[2]))
quantile_cover$wr_sd_1 = quantile_coverage(scalars_out$wr_sd_1,param_settings_in$wr_sd[2])

print("wr SD, dual")
print(quantile(scalars_out$wr_sd_2,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$wr_sd_2,param_settings_in$wr_sd[3]))
quantile_cover$wr_sd_2 = quantile_coverage(scalars_out$wr_sd_2,param_settings_in$wr_sd[3])

print("sigma")
print(quantile(scalars_out$sigma,probs=c(0.025,0.5,0.975)))
print(quantile_coverage(scalars_out$sigma,param_settings_in$sigma))
quantile_cover$sigma = quantile_coverage(scalars_out$sigma,param_settings_in$sigma)

acp_pr <- function(mcmc_vec){
  return(mean(mcmc_vec[1:(length(mcmc_vec)-1)] != mcmc_vec[2:(length(mcmc_vec))],na.rm=T))
}

scalar_acp_pr <- apply(scalars_out,2,acp_pr)

print("Acceptance Probabilities")
print(scalar_acp_pr)

print("ESS")
print(apply(scalars_out,2,effectiveSize))

sink(file=NULL)

save(quantile_cover, file = "quantile_coverage.RData")

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

model_estimate <- rep(NA,ncol(model_out))

pdf(file="Model_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))

for(i in 1:ncol(model_out)){
  model_p <- round(table(c(model_out[,i],1,2,3))/(length(model_out[,i])+3),3)
  
  plot(model_out[,i],main=paste0("Model, ",index_id[i],", p=(",model_p[1],",",model_p[2],",",model_p[3],")"),type="l",xlab="Iteration",ylab="Model")
  abline(h=individual_data$model_true[i],col="blue")
  
  model_estimate[i] <- which.max(model_p)
}

par(mfrow=c(1,1))

dev.off()

true_model <- individual_data$model_true

model_correct_summary <- table(true_model,model_estimate)
p_model_correct <- sum(diag(model_correct_summary))/sum(model_correct_summary)

save(model_estimate,model_correct_summary,p_model_correct,file="model_correct_summary.RData")

print("WP Traceplots")

pdf(file="WP_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))
for(i in 1:ncol(wp_out)){
  plot(wp_out[,i],main=paste0("wp, individual ",index_id[i]),type="l",xlab="Iteration",ylab="wp")
  abline(h=individual_data$wp_true[i])
}
par(mfrow=c(1,1))

dev.off()

print("TP Traceplots")

pdf(file="TP_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))
for(i in 1:ncol(tp_out)){
  plot(tp_out[,i],main=paste0("Tp, individual ",index_id[i]),type="l",xlab="Iteration",ylab="tp")
  abline(h=individual_data$tp_true[i])
}
par(mfrow=c(1,1))

dev.off()

print("DP Traceplots")

pdf(file="DP_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))
for(i in 1:ncol(dp_out)){
  plot(dp_out[,i],main=paste0("dp, individual ",index_id[i]),type="l",xlab="Iteration",ylab="dp")
  abline(h=individual_data$dp_true[i])
}
par(mfrow=c(1,1))

dev.off()

print("WR Traceplots")

pdf(file="WR_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))
for(i in 1:ncol(wr_out)){
  plot(wr_out[,i],main=paste0("wr, individual ",index_id[i]),type="l",xlab="Iteration",ylab="wr")
  abline(h=individual_data$wr_true[i])
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
  
  model_infer <- which.max(table(c(model_out[,i],1,2,3)))
  
  model_true <- individual_data$model_true[i]
  
  line_col <- ifelse(model_infer == model_true,"black","red4")
  
  plot_col <- c("blue","red","purple")[individual_data$subtype[i]+1]
  
  # plot_col <- "blue"
  
  select_model_iter <- model_out[,i] == model_infer
  
  wp_quantile <- quantile(wp_out[select_model_iter,i],probs = c(0.025,0.5,0.975))
  tp_quantile <- quantile(tp_out[select_model_iter,i],probs = c(0.025,0.5,0.975))
  dp_quantile <- quantile(dp_out[select_model_iter,i],probs = c(0.025,0.5,0.975))
  wr_quantile <- quantile(wr_out[select_model_iter,i],probs = c(0.025,0.5,0.975))
  
  tp_med <- median(tp_out[select_model_iter,i])
  
  plot(plot_dat$time,plot_dat$viral_load,pch=19,main=paste0("Observed Data, ID ",i-1),col=plot_col,
       xlim=c(min(-wp_quantile[3],plot_dat$time),max(wr_quantile[3],plot_dat$time)),
       ylim=c(0,max(dp_quantile[3],plot_dat$time)))
  
  if(model_infer != 3){
    lines(c(tp_med-wp_quantile[2],tp_med),c(0,dp_quantile[2]),col=line_col)
    lines(c(tp_med-wp_quantile[1],tp_med),c(0,dp_quantile[1]),lty="dashed",col=line_col)
    lines(c(tp_med-wp_quantile[3],tp_med),c(0,dp_quantile[3]),lty="dashed",col=line_col)
  }
  if(model_infer != 1){
    lines(c(tp_med,tp_med+wr_quantile[2]),c(dp_quantile[2],0),col=line_col)
    lines(c(tp_med,tp_med+wr_quantile[1]),c(dp_quantile[1],0),lty="dashed",col=line_col)
    lines(c(tp_med,tp_med+wr_quantile[3]),c(dp_quantile[3],0),lty="dashed",col=line_col)
  }
  
  
}

par(mfrow=c(1,1))
dev.off()


### Violin plots
### Parameter violin plots
lod <- settings$lod

shared_sd = T

wp_mean_violin_data <- data.frame(subtype=factor(c(rep("Subtype 1",nrow(scalars_out)),rep("Subtype 2",nrow(scalars_out)),rep("Subtype 3",nrow(scalars_out))),levels=c("Subtype 1","Subtype 2","Subtype 3")),
                                  value=c(scalars_out$wp_mean_0,scalars_out$wp_mean_1,scalars_out$wp_mean_2))

wp_mean_violin <- ggplot(data=wp_mean_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=2.5) +
  scale_fill_manual(values=c("skyblue2","coral","purple3"),guide="none") +
  labs(x="",y="Mean (Days)",fill="Subtype") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
  annotate("text", x = 1, y = param_settings_in$wp_mean[1], label = "*", size = 8, color = "black", fontface = "bold") +
  annotate("text", x = 2, y = param_settings_in$wp_mean[1], label = "*", size = 8, color = "black", fontface = "bold") +
  annotate("text", x = 3, y = param_settings_in$wp_mean[1], label = "*", size = 8, color = "black", fontface = "bold") +
  ggtitle("Mean Proliferation Duration")

wp_mean_violin

dp_mean_violin_data <- data.frame(subtype=factor(c(rep("Subtype 1",nrow(scalars_out)),rep("Subtype 2",nrow(scalars_out)),rep("Subtype 3",nrow(scalars_out))),levels=c("Subtype 1","Subtype 2","Subtype 3")),
                                  value=c(scalars_out$dp_mean_0,scalars_out$dp_mean_1,scalars_out$dp_mean_2))

dp_mean_violin <- ggplot(data=dp_mean_violin_data,aes(x=factor(subtype),y=lod-value,fill=subtype,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=4) +
  scale_fill_manual(values=c("skyblue2","coral","purple3"),guide="none") +
  scale_y_reverse() +
  labs(x="",y="Mean Peak Ct",fill="Subtype") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
  annotate("text", x = 1, y = lod - param_settings_in$dp_mean[1], label = "*", size = 8, color = "black", fontface = "bold") +
  annotate("text", x = 2, y = lod - param_settings_in$dp_mean[1], label = "*", size = 8, color = "black", fontface = "bold") +
  annotate("text", x = 3, y = lod - param_settings_in$dp_mean[1], label = "*", size = 8, color = "black", fontface = "bold") +
  ggtitle("Mean Viral Load Peak")

dp_mean_violin

wr_mean_violin_data <- data.frame(subtype=factor(c(rep("Subtype 1",nrow(scalars_out)),rep("Subtype 2",nrow(scalars_out)),rep("Subtype 3",nrow(scalars_out))),levels=c("Subtype 1","Subtype 2","Subtype 3")),
                                  value=c(scalars_out$wr_mean_0,scalars_out$wr_mean_1,scalars_out$wr_mean_2))

wr_mean_violin <- ggplot(data=wr_mean_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=4) +
  scale_fill_manual(values=c("skyblue2","coral","purple3"),guide="none") +
  labs(x="",y="Mean (Days)",fill="Subtype") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
  annotate("text", x = 1, y = param_settings_in$wr_mean[1], label = "*", size = 8, color = "black", fontface = "bold") +
  annotate("text", x = 2, y = param_settings_in$wr_mean[2], label = "*", size = 8, color = "black", fontface = "bold") +
  annotate("text", x = 3, y = param_settings_in$wr_mean[3], label = "*", size = 8, color = "black", fontface = "bold") +
  ggtitle("Mean Clearance Duration")

wr_mean_violin

if(shared_sd){
  wp_sd_violin_data <- data.frame(subtype=factor(rep("All Subtypes",nrow(scalars_out)),levels=c("All Subtypes")),
                                  value=c(scalars_out$wp_sd_0))
  
  wp_sd_violin <- ggplot(data=wp_sd_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
    geom_violin(trim=F,linewidth=0.2,adjust=2.5,width=0.35) +
    scale_fill_manual(values=c("grey30"),guide="none") +
    labs(x="",y="SD Proliferation",fill="Subtype") +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size=11)) +
    annotate("text", x = 1, y = param_settings_in$wp_sd[1], label = "*", size = 8, color = "black", fontface = "bold") +
    ggtitle("SD, Proliferation Duration")
} else{
  wp_sd_violin_data <- data.frame(subtype=factor(c(rep("Subtype 1",nrow(scalars_out)),rep("Subtype 2",nrow(scalars_out)),rep("Subtype 3",nrow(scalars_out))),levels=c("Subtype 1","Subtype 2","Subtype 3")),
                                  value=c(scalars_out$wp_sd_0,scalars_out$wp_sd_1,scalars_out$wp_sd_2))
  
  wp_sd_violin <- ggplot(data=wp_sd_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
    geom_violin(trim=F,linewidth=0.2,adjust=2.5) +
    scale_fill_manual(values=c("skyblue2","coral","purple3"),guide="none") +
    labs(x="",y="SD Proliferation",fill="Subtype") +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size=11)) +
    annotate("text", x = 1, y = param_settings_in$wp_sd[1], label = "*", size = 8, color = "black", fontface = "bold") +
    ggtitle("SD, Proliferation Duration")
}

wp_sd_violin

if(shared_sd){
  dp_sd_violin_data <- data.frame(subtype=factor(rep("All Subtypes",nrow(scalars_out)),levels=c("All Subtypes")),
                                  value=c(scalars_out$dp_sd_0))
  
  dp_sd_violin <- ggplot(data=dp_sd_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
    geom_violin(trim=F,linewidth=0.2,adjust=2.5,width=0.35) +
    scale_fill_manual(values=c("grey30"),guide="none") +
    labs(x="",y="SD Peak",fill="Subtype") +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size=11)) +
    annotate("text", x = 1, y = param_settings_in$dp_sd[1], label = "*", size = 8, color = "black", fontface = "bold") +
    ggtitle("SD, Viral Load Peak")
} else{
  dp_sd_violin_data <- data.frame(subtype=factor(c(rep("Subtype 1",nrow(scalars_out)),rep("Subtype 2",nrow(scalars_out)),rep("Subtype 3",nrow(scalars_out))),levels=c("Subtype 1","Subtype 2","Subtype 3")),
                                  value=c(scalars_out$dp_sd_0,scalars_out$dp_sd_1,scalars_out$dp_sd_2))
  
  dp_sd_violin <- ggplot(data=dp_sd_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
    geom_violin(trim=F,linewidth=0.2,adjust=2.5) +
    scale_fill_manual(values=c("skyblue2","coral","purple3"),guide="none") +
    labs(x="",y="SD Peak",fill="Subtype") +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size=11)) +
    ggtitle("SD, Viral Load Peak")
}

dp_sd_violin

if(shared_sd){
  wr_sd_violin_data <- data.frame(subtype=factor(rep("All Subtypes",nrow(scalars_out)),levels=c("All Subtypes")),
                                  value=c(scalars_out$wr_sd_0))
  
  wr_sd_violin <- ggplot(data=wr_sd_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
    geom_violin(trim=F,linewidth=0.2,adjust=2.5,width=0.35) +
    scale_fill_manual(values=c("grey30"),guide="none") +
    labs(x="",y="SD Clearance",fill="Subtype") +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size=11)) +
    annotate("text", x = 1, y = param_settings_in$wr_sd[1], label = "*", size = 8, color = "black", fontface = "bold") +
    ggtitle("SD, Clearance Duration")
} else{
  wr_sd_violin_data <- data.frame(subtype=factor(c(rep("Subtype 1",nrow(scalars_out)),rep("Subtype 2",nrow(scalars_out)),rep("Subtype 3",nrow(scalars_out))),levels=c("Subtype 1","Subtype 2","Subtype 3")),
                                  value=c(scalars_out$wr_sd_0,scalars_out$wr_sd_1,scalars_out$wr_sd_2))
  
  wr_sd_violin <- ggplot(data=wr_sd_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
    geom_violin(trim=F,linewidth=0.2,adjust=2.5) +
    scale_fill_manual(values=c("skyblue2","coral","purple3"),guide="none") +
    labs(x="",y="SD Clearance",fill="Subtype") +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size=11)) +
    ggtitle("SD, Clearance Duration")
}


wr_sd_violin

sigma_sd_violin_data <- data.frame(value = scalars_out$sigma,
                                    subtype = factor("All Subtypes"))

sigma_violin <- ggplot(data=sigma_sd_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=2.5,width=0.35) +
  scale_fill_manual(values=c("grey30"),guide="none") +
  labs(x="",y="SD",fill="Subtype") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
  annotate("text", x = 1, y = param_settings_in$sigma, label = "*", size = 8, color = "black", fontface = "bold") +
  ggtitle("SD, Viral Observations")

sigma_violin

posterior_violin_plots_trajectories <- ggarrange(wp_mean_violin,
                                                 dp_mean_violin,
                                                 wr_mean_violin,
                                                 wp_sd_violin,
                                                 dp_sd_violin,
                                                 wr_sd_violin,
                                                 NULL,
                                                 sigma_violin,
                                                 NULL,
                                                 nrow=3,ncol=3,
                                                 labels = c("(a)","(b)","(c)","(d)","(e)","(f)","","(g)",""))

annotate_figure(posterior_violin_plots_trajectories, top = text_grob("Posterior Distributions", 
                                                                     face = "bold", size = 14))

png(filename="posterior_violins_mean_sd_multichain_sim.png",width=2809,height=3350,res=300)
print(posterior_violin_plots_trajectories)
dev.off()

lod <- settings$lod

t_min <- -7
t_max <- 20

plot_dat_H1N1 <- data.frame(t_val = c(t_min,
                                      -quantile(scalars_out$wp_mean_0,probs=0.975),
                                      -quantile(scalars_out$wp_mean_0,probs=0.5),
                                      -quantile(scalars_out$wp_mean_0,probs=0.025),
                                      0,
                                      quantile(scalars_out$wr_mean_0,probs=0.025),
                                      quantile(scalars_out$wr_mean_0,probs=0.5),
                                      quantile(scalars_out$wr_mean_0,probs=0.975),
                                      t_max),
                            med_val = c(0,
                                        0,
                                        0,
                                        mu(-quantile(scalars_out$wp_mean_0,probs=0.025),
                                           quantile(scalars_out$wp_mean_0,probs=0.5),
                                           0,
                                           quantile(scalars_out$dp_mean_0,probs=0.5),
                                           quantile(scalars_out$wr_mean_0,probs=0.5)),
                                        quantile(scalars_out$dp_mean_0,probs=0.5),
                                        mu(quantile(scalars_out$wr_mean_0,probs=0.025),
                                           quantile(scalars_out$wp_mean_0,probs=0.5),
                                           0,
                                           quantile(scalars_out$dp_mean_0,probs=0.5),
                                           quantile(scalars_out$wr_mean_0,probs=0.5)),
                                        0,
                                        0,
                                        0),
                            l_val = c(0,
                                      0,
                                      0,
                                      0,
                                      quantile(scalars_out$dp_mean_0,probs=0.025),
                                      0,
                                      0,
                                      0,
                                      0),
                            h_val = c(0,
                                      0,
                                      mu(-quantile(scalars_out$wp_mean_0,probs=0.5),
                                         quantile(scalars_out$wp_mean_0,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_0,probs=0.975),
                                         quantile(scalars_out$wr_mean_0,probs=0.975)),
                                      mu(-quantile(scalars_out$wp_mean_0,probs=0.025),
                                         quantile(scalars_out$wp_mean_0,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_0,probs=0.975),
                                         quantile(scalars_out$wr_mean_0,probs=0.975)),
                                      quantile(scalars_out$dp_mean_0,probs=0.975),
                                      mu(quantile(scalars_out$wr_mean_0,probs=0.025),
                                         quantile(scalars_out$wp_mean_0,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_0,probs=0.975),
                                         quantile(scalars_out$wr_mean_0,probs=0.975)),
                                      mu(quantile(scalars_out$wr_mean_0,probs=0.5),
                                         quantile(scalars_out$wp_mean_0,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_0,probs=0.975),
                                         quantile(scalars_out$wr_mean_0,probs=0.975)),
                                      0,
                                      0))

plot_dat_H3N2 <- data.frame(t_val = c(t_min,
                                      -quantile(scalars_out$wp_mean_1,probs=0.975),
                                      -quantile(scalars_out$wp_mean_1,probs=0.5),
                                      -quantile(scalars_out$wp_mean_1,probs=0.025),
                                      0,
                                      quantile(scalars_out$wr_mean_1,probs=0.025),
                                      quantile(scalars_out$wr_mean_1,probs=0.5),
                                      quantile(scalars_out$wr_mean_1,probs=0.975),
                                      t_max),
                            med_val = c(0,
                                        0,
                                        0,
                                        mu(-quantile(scalars_out$wp_mean_1,probs=0.025),
                                           quantile(scalars_out$wp_mean_1,probs=0.5),
                                           0,
                                           quantile(scalars_out$dp_mean_1,probs=0.5),
                                           quantile(scalars_out$wr_mean_1,probs=0.5)),
                                        quantile(scalars_out$dp_mean_1,probs=0.5),
                                        mu(quantile(scalars_out$wr_mean_1,probs=0.025),
                                           quantile(scalars_out$wp_mean_1,probs=0.5),
                                           0,
                                           quantile(scalars_out$dp_mean_1,probs=0.5),
                                           quantile(scalars_out$wr_mean_1,probs=0.5)),
                                        0,
                                        0,
                                        0),
                            l_val = c(0,
                                      0,
                                      0,
                                      0,
                                      quantile(scalars_out$dp_mean_1,probs=0.025),
                                      0,
                                      0,
                                      0,
                                      0),
                            h_val = c(0,
                                      0,
                                      mu(-quantile(scalars_out$wp_mean_1,probs=0.5),
                                         quantile(scalars_out$wp_mean_1,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_1,probs=0.975),
                                         quantile(scalars_out$wr_mean_1,probs=0.975)),
                                      mu(-quantile(scalars_out$wp_mean_1,probs=0.025),
                                         quantile(scalars_out$wp_mean_1,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_1,probs=0.975),
                                         quantile(scalars_out$wr_mean_1,probs=0.975)),
                                      quantile(scalars_out$dp_mean_1,probs=0.975),
                                      mu(quantile(scalars_out$wr_mean_1,probs=0.025),
                                         quantile(scalars_out$wp_mean_1,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_1,probs=0.975),
                                         quantile(scalars_out$wr_mean_1,probs=0.975)),
                                      mu(quantile(scalars_out$wr_mean_1,probs=0.5),
                                         quantile(scalars_out$wp_mean_1,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_1,probs=0.975),
                                         quantile(scalars_out$wr_mean_1,probs=0.975)),
                                      0,
                                      0))

plot_dat_dual <- data.frame(t_val = c(t_min,
                                      -quantile(scalars_out$wp_mean_2,probs=0.975),
                                      -quantile(scalars_out$wp_mean_2,probs=0.5),
                                      -quantile(scalars_out$wp_mean_2,probs=0.025),
                                      0,
                                      quantile(scalars_out$wr_mean_2,probs=0.025),
                                      quantile(scalars_out$wr_mean_2,probs=0.5),
                                      quantile(scalars_out$wr_mean_2,probs=0.975),
                                      t_max),
                            med_val = c(0,
                                        0,
                                        0,
                                        mu(-quantile(scalars_out$wp_mean_2,probs=0.025),
                                           quantile(scalars_out$wp_mean_2,probs=0.5),
                                           0,
                                           quantile(scalars_out$dp_mean_2,probs=0.5),
                                           quantile(scalars_out$wr_mean_2,probs=0.5)),
                                        quantile(scalars_out$dp_mean_2,probs=0.5),
                                        mu(quantile(scalars_out$wr_mean_2,probs=0.025),
                                           quantile(scalars_out$wp_mean_2,probs=0.5),
                                           0,
                                           quantile(scalars_out$dp_mean_2,probs=0.5),
                                           quantile(scalars_out$wr_mean_2,probs=0.5)),
                                        0,
                                        0,
                                        0),
                            l_val = c(0,
                                      0,
                                      0,
                                      0,
                                      quantile(scalars_out$dp_mean_2,probs=0.025),
                                      0,
                                      0,
                                      0,
                                      0),
                            h_val = c(0,
                                      0,
                                      mu(-quantile(scalars_out$wp_mean_2,probs=0.5),
                                         quantile(scalars_out$wp_mean_2,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_2,probs=0.975),
                                         quantile(scalars_out$wr_mean_2,probs=0.975)),
                                      mu(-quantile(scalars_out$wp_mean_2,probs=0.025),
                                         quantile(scalars_out$wp_mean_2,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_2,probs=0.975),
                                         quantile(scalars_out$wr_mean_2,probs=0.975)),
                                      quantile(scalars_out$dp_mean_2,probs=0.975),
                                      mu(quantile(scalars_out$wr_mean_2,probs=0.025),
                                         quantile(scalars_out$wp_mean_2,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_2,probs=0.975),
                                         quantile(scalars_out$wr_mean_2,probs=0.975)),
                                      mu(quantile(scalars_out$wr_mean_2,probs=0.5),
                                         quantile(scalars_out$wp_mean_2,probs=0.975),
                                         0,
                                         quantile(scalars_out$dp_mean_2,probs=0.975),
                                         quantile(scalars_out$wr_mean_2,probs=0.975)),
                                      0,
                                      0))

plot_y_min <- lod
plot_y_max <- min(lod-c(plot_dat_H1N1$h_val,plot_dat_H3N2$h_val,plot_dat_dual$h_val))

H1N1_trajectory_plot <- plot_dat_H1N1 %>% ggplot(aes(x=t_val,y=lod-med_val)) +
  geom_ribbon(aes(ymin=lod-l_val,ymax=lod-h_val), alpha=0.5, linewidth = 0, fill="skyblue2") +
  geom_line(col="skyblue2") +
  labs(title = "Mean Viral Trajectory, H1N1",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Ct") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

print(H1N1_trajectory_plot)

H3N2_trajectory_plot <- plot_dat_H3N2 %>% ggplot(aes(x=t_val,y=lod-med_val)) +
  geom_ribbon(aes(ymin=lod-l_val,ymax=lod-h_val), alpha=0.5, linewidth = 0, fill="coral") +
  geom_line(col="coral") +
  labs(title = "Mean Viral Trajectory, H3N2",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Ct") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

print(H3N2_trajectory_plot)

dual_trajectory_plot <- plot_dat_dual %>% ggplot(aes(x=t_val,y=lod-med_val)) +
  geom_ribbon(aes(ymin=lod-l_val,ymax=lod-h_val), alpha=0.5, linewidth = 0, fill="purple3") +
  geom_line(col="purple3") +
  labs(title = "Mean Viral Trajectory, Dual",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Ct") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

print(dual_trajectory_plot)

n_iter_total <- nrow(scalars_out)
# end_iter <- (1:(length(seeds)))/length(seeds)*n_iter_total
# begin_iter <- c(1,end_iter[1:(length(seeds)-1)]-1)
# 
# mcmc_multi_list <- as.mcmc.list(as.mcmc(scalars_out[begin_iter[1]:end_iter[1],c("wp_mean_0","wp_mean_1","wp_mean_2","wp_sd_0","wp_sd_1","wp_sd_2","dp_mean_0","dp_mean_1","dp_mean_2","dp_sd_0","dp_sd_1","dp_sd_2","wr_mean_0","wr_mean_1","wr_mean_2","wr_sd_0","wr_sd_1","wr_sd_2","sigma")]))
# 
# for(i in 2:length(seeds)){
#   mcmc_multi_list[[i]] <- as.mcmc(scalars_out[begin_iter[i]:end_iter[i],c("wp_mean_0","wp_mean_1","wp_mean_2","wp_sd_0","wp_sd_1","wp_sd_2","dp_mean_0","dp_mean_1","dp_mean_2","dp_sd_0","dp_sd_1","dp_sd_2","wr_mean_0","wr_mean_1","wr_mean_2","wr_sd_0","wr_sd_1","wr_sd_2","sigma")])
# }
# 
# sink(file="gelman.txt")
# print(gelman.diag(mcmc_multi_list,multivariate = F))
# sink(file=NULL)


### Sample mu from posteriors, make data plot
xmin <- -5
xmax <- 20

xval <- seq(from=xmin,to=xmax,length.out=200)
n_sample <- 10000

ydat_h1n1 <- matrix(data=NA,ncol=length(xval),nrow=n_sample)
trajectory_dat_h1n1 <- matrix(data=NA,ncol=3,nrow=n_sample)

# trajectory_dat_h1n1[,1] <- rnorm(n_sample,median(scalars_out$wp_mean_0),median(scalars_out$wp_sd_0))
# trajectory_dat_h1n1[,2] <- rnorm(n_sample,median(scalars_out$dp_mean_0),median(scalars_out$dp_sd_0))
# trajectory_dat_h1n1[,3] <- rnorm(n_sample,median(scalars_out$wr_mean_0),median(scalars_out$wr_sd_0))

trajectory_dat_h1n1[,1] <- sample(scalars_out$wp_mean_0,n_sample,replace=T)
trajectory_dat_h1n1[,2] <- sample(scalars_out$dp_mean_0,n_sample,replace=T)
trajectory_dat_h1n1[,3] <- sample(scalars_out$wr_mean_0,n_sample,replace=T)

for(i in 1:nrow(trajectory_dat_h1n1)){
  for(j in 1:length(xval)){
    ydat_h1n1[i,j] <- ViralLoadRJMCMC::mu(xval[j],trajectory_dat_h1n1[i,1],0,trajectory_dat_h1n1[i,2],trajectory_dat_h1n1[i,3])
  }
}

sample_trajectory_h1n1 <- apply(ydat_h1n1,2,quantile,probs=c(0.025,0.5,0.975))

trajectory_sample_dat_h1n1 <- data.frame(timeval=xval,
                                         median=sample_trajectory_h1n1[2,],
                                         low=sample_trajectory_h1n1[1,],
                                         hi=sample_trajectory_h1n1[3,])

plot_y_min <- settings$lod
plot_y_max <- settings$lod - max(trajectory_sample_dat_h1n1$hi)

### H3N2
ydat_h3n2 <- matrix(data=NA,ncol=length(xval),nrow=n_sample)
trajectory_dat_h3n2 <- matrix(data=NA,ncol=3,nrow=n_sample)

# trajectory_dat_h3n2[,1] <- rnorm(n_sample,median(scalars_out$wp_mean_1),median(scalars_out$wp_sd_1))
# trajectory_dat_h3n2[,2] <- rnorm(n_sample,median(scalars_out$dp_mean_1),median(scalars_out$dp_sd_1))
# trajectory_dat_h3n2[,3] <- rnorm(n_sample,median(scalars_out$wr_mean_1),median(scalars_out$wr_sd_1))

trajectory_dat_h3n2[,1] <- sample(scalars_out$wp_mean_1,n_sample,replace=T)
trajectory_dat_h3n2[,2] <- sample(scalars_out$dp_mean_1,n_sample,replace=T)
trajectory_dat_h3n2[,3] <- sample(scalars_out$wr_mean_1,n_sample,replace=T)

for(i in 1:nrow(trajectory_dat_h3n2)){
  for(j in 1:length(xval)){
    ydat_h3n2[i,j] <- ViralLoadRJMCMC::mu(xval[j],trajectory_dat_h3n2[i,1],0,trajectory_dat_h3n2[i,2],trajectory_dat_h3n2[i,3])
  }
}

sample_trajectory_h3n2 <- apply(ydat_h3n2,2,quantile,probs=c(0.025,0.5,0.975))

trajectory_sample_dat_h3n2 <- data.frame(timeval=xval,
                                         median=sample_trajectory_h3n2[2,],
                                         low=sample_trajectory_h3n2[1,],
                                         hi=sample_trajectory_h3n2[3,])

plot_y_min <- settings$lod
plot_y_max <- settings$lod - max(trajectory_sample_dat_h3n2$hi)

### Dual
ydat_dual <- matrix(data=NA,ncol=length(xval),nrow=n_sample)
trajectory_dat_dual <- matrix(data=NA,ncol=3,nrow=n_sample)

# trajectory_dat_dual[,1] <- rnorm(n_sample,median(scalars_out$wp_mean_2),median(scalars_out$wp_sd_2))
# trajectory_dat_dual[,2] <- rnorm(n_sample,median(scalars_out$dp_mean_2),median(scalars_out$dp_sd_2))
# trajectory_dat_dual[,3] <- rnorm(n_sample,median(scalars_out$wr_mean_2),median(scalars_out$wr_sd_2))

trajectory_dat_dual[,1] <- sample(scalars_out$wp_mean_2,n_sample,replace=T)
trajectory_dat_dual[,2] <- sample(scalars_out$dp_mean_2,n_sample,replace=T)
trajectory_dat_dual[,3] <- sample(scalars_out$wr_mean_2,n_sample,replace=T)

for(i in 1:nrow(trajectory_dat_dual)){
  for(j in 1:length(xval)){
    ydat_dual[i,j] <- ViralLoadRJMCMC::mu(xval[j],trajectory_dat_dual[i,1],0,trajectory_dat_dual[i,2],trajectory_dat_dual[i,3])
  }
}

sample_trajectory_dual <- apply(ydat_dual,2,quantile,probs=c(0.025,0.5,0.975))

trajectory_sample_dat_dual <- data.frame(timeval=xval,
                                         median=sample_trajectory_dual[2,],
                                         low=sample_trajectory_dual[1,],
                                         hi=sample_trajectory_dual[3,])

plot_y_min <- settings$lod
plot_y_max <- settings$lod - max(trajectory_sample_dat_h1n1$hi,trajectory_sample_dat_h3n2$hi,trajectory_sample_dat_dual$hi)

trajectory_sample_plot_h1n1 <- trajectory_sample_dat_h1n1 %>% ggplot(aes(x=timeval,y=settings$lod-median)) +
  geom_ribbon(aes(ymin=settings$lod-low,ymax=settings$lod-hi), alpha=0.5, linewidth = 0, fill="skyblue2") +
  geom_line(col="skyblue2") +
  labs(title = "Posterior Sampled Trajectories, H1N1",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Ct") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

trajectory_sample_plot_h3n2 <- trajectory_sample_dat_h3n2 %>% ggplot(aes(x=timeval,y=settings$lod-median)) +
  geom_ribbon(aes(ymin=settings$lod-low,ymax=settings$lod-hi), alpha=0.5, linewidth = 0, fill="coral") +
  geom_line(col="coral") +
  labs(title = "Posterior Sampled Trajectories, H3N2",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Ct") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

trajectory_sample_plot_dual <- trajectory_sample_dat_dual %>% ggplot(aes(x=timeval,y=settings$lod-median)) +
  geom_ribbon(aes(ymin=settings$lod-low,ymax=settings$lod-hi), alpha=0.5, linewidth = 0, fill="purple3") +
  geom_line(col="purple3") +
  labs(title = "Posterior Sampled Trajectories, Dual",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Ct") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

png(file="sample_trajectory_flu.png",width=1100,height=1300,res=200)

print(ggarrange(trajectory_sample_plot_h1n1,trajectory_sample_plot_h3n2,trajectory_sample_plot_dual,ncol=1))

dev.off()

### Summarize true wp values etc
wp_interval <- apply(wp_out,2,quantile,probs=c(0.025,0.5,0.975))
dp_interval <- apply(dp_out,2,quantile,probs=c(0.025,0.5,0.975))
tp_interval <- apply(tp_out,2,quantile,probs=c(0.025,0.5,0.975))
wr_interval <- apply(wr_out,2,quantile,probs=c(0.025,0.5,0.975))

# hdi_median <- function(x, cint=0.95){
#   hdi_x <- bayestestR::hdi(x, ci = cint)
#   med_x <- median(x)
#   
#   return(unlist(c(hdi_x[2], med_x, hdi_x[3])))
# }
# 
# wp_interval <- apply(wp_out,2,hdi_median)
# dp_interval <- apply(dp_out,2,hdi_median)
# tp_interval <- apply(tp_out,2,hdi_median)
# wr_interval <- apply(wr_out,2,hdi_median)

who_wp <- which(individual_data$model_true %in% c(1,2))
mean(wp_interval[1,who_wp] <= individual_data$wp_true[who_wp] & wp_interval[3,who_wp] >= individual_data$wp_true[who_wp])
mean(dp_interval[1,] <= individual_data$dp_true & dp_interval[3,] >= individual_data$dp_true)
mean(tp_interval[1,] <= individual_data$tp_true & tp_interval[3,] >= individual_data$tp_true)
who_wr <- which(individual_data$model_true %in% c(2,3))
mean(wr_interval[1,who_wr] <= individual_data$wr_true[who_wr] & wr_interval[3,who_wr] >= individual_data$wr_true[who_wr])

wp_intervals <- matrix(data=NA, nrow = 2, ncol = ncol(wp_out))
dp_intervals <- matrix(data=NA, nrow = 2, ncol = ncol(dp_out))
tp_intervals <- matrix(data=NA, nrow = 2, ncol = ncol(tp_out))
wr_intervals <- matrix(data=NA, nrow = 2, ncol = ncol(wr_out))

wp_interval_cover <- rep(NA, ncol(wp_out))
dp_interval_cover <- rep(NA, ncol(dp_out))
tp_interval_cover <- rep(NA, ncol(tp_out))
wr_interval_cover <- rep(NA, ncol(wr_out))

for(i in 1:ncol(wp_out)){
  model_i <- individual_data$model_true[i]
  wp_i <- individual_data$wp_true[i]
  dp_i <- individual_data$dp_true[i]
  tp_i <- individual_data$tp_true[i]
  wr_i <- individual_data$wr_true[i]
  
  true_model_iter <- which(model_out[,i] == model_i)
  print(mean(model_out[,i] == model_i))
  if(length(true_model_iter)==0){
    next
  }
  
  if(model_i %in% c(1)){
    wp_interval <- quantile(wp_out[true_model_iter,i], probs=c(0.025, 0.975))
    
    wp_intervals[,i] <- wp_interval
    
    wp_interval_cover[i] <- wp_interval[1] <= wp_i & wp_interval[2] >= wp_i
  }
  
  if(model_i %in% c(3)){
    wr_interval <- quantile(wr_out[true_model_iter,i], probs=c(0.025, 0.975))
    
    wr_intervals[,i] <- wr_interval
    
    wr_interval_cover[i] <- wr_interval[1] <= wr_i & wr_interval[2] >= wr_i
  }
  
  dp_interval <- quantile(dp_out[true_model_iter,i], probs=c(0.025, 0.975))
  
  dp_intervals[,i] <- dp_interval
  
  dp_interval_cover[i] <- dp_interval[1] <= dp_i & dp_interval[2] >= dp_i
  
  tp_interval <- quantile(tp_out[true_model_iter,i], probs=c(0.025, 0.975))
  
  tp_intervals[,i] <- tp_interval
  
  tp_interval_cover[i] <- tp_interval[1] <= tp_i & tp_interval[2] >= tp_i
}

mean(wp_interval_cover, na.rm=T)
mean(dp_interval_cover, na.rm=T)
mean(tp_interval_cover, na.rm=T)
mean(wr_interval_cover, na.rm=T)

