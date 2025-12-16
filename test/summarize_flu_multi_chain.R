library(Rcpp)
library(dplyr)
library(truncnorm)
library(ggpubr)
library(ggplot2)
library(coda)
library(ViralLoadRJMCMC)

dec.precision <- 6

pct_burnin_begin <- 0.5
pct_burnin_end <- 0.99

results_dir <- "~/Documents/Research/Within-Host/RJMCMC_Results"

setwd(results_dir)

seeds <- c(1111:1115)
# seeds <- c(1111:1114)

shared_sd <- T

# seeds <- c(1111,1112,1114,1115)

analysis_dirs <- paste0("./analysis_flu_seed_",seeds,"/")

setwd(analysis_dirs[1])

load("individual_data_in.RData")
load("priors.RData")
load("settings_in.RData")
load("viral_data_in.RData")

index_id <- 0:(nrow(individual_data)-1)

scalars_out_raw <- read.csv("./output/scalars_out.csv",header=T)
n_burnin_begin <- round(nrow(scalars_out_raw)*pct_burnin_begin)
n_burnin_end <- round(nrow(scalars_out_raw)*pct_burnin_end)

scalars_out <- scalars_out_raw[n_burnin_begin:n_burnin_end,]
rm(scalars_out_raw)

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

setwd(results_dir)

for(i in 2:length(analysis_dirs)){
  setwd(analysis_dirs[i])
  
  scalars_out_raw <- read.csv("./output/scalars_out.csv",header=T)
  n_burnin_begin <- round(nrow(scalars_out_raw)*pct_burnin_begin)
  n_burnin_end <- round(nrow(scalars_out_raw)*pct_burnin_end)
  
  scalars_out <- rbind(scalars_out,scalars_out_raw[n_burnin_begin:n_burnin_end,])
  rm(scalars_out_raw)
  
  model_out_raw <- read.csv("./output/model_out.csv",header=F)
  n_burnin_begin <- round(nrow(model_out_raw)*pct_burnin_begin)
  n_burnin_end <- round(nrow(model_out_raw)*pct_burnin_end)
  
  model_out <- rbind(model_out,model_out_raw[n_burnin_begin:n_burnin_end,])
  rm(model_out_raw)
  
  wp_out_raw <- read.csv("./output/wp_out.csv",header=F)
  n_burnin_begin <- round(nrow(wp_out_raw)*pct_burnin_begin)
  n_burnin_end <- round(nrow(wp_out_raw)*pct_burnin_end)
  
  wp_out <- rbind(wp_out,wp_out_raw[n_burnin_begin:n_burnin_end,])
  rm(wp_out_raw)
  
  tp_out_raw <- read.csv("./output/tp_out.csv",header=F)
  n_burnin_begin <- round(nrow(tp_out_raw)*pct_burnin_begin)
  n_burnin_end <- round(nrow(tp_out_raw)*pct_burnin_end)
  
  tp_out <- rbind(tp_out,tp_out_raw[n_burnin_begin:n_burnin_end,])
  rm(tp_out_raw)
  
  dp_out_raw <- read.csv("./output/dp_out.csv",header=F)
  n_burnin_begin <- round(nrow(dp_out_raw)*pct_burnin_begin)
  n_burnin_end <- round(nrow(dp_out_raw)*pct_burnin_end)
  
  dp_out <- rbind(dp_out,dp_out_raw[n_burnin_begin:n_burnin_end,])
  rm(dp_out_raw)
  
  wr_out_raw <- read.csv("./output/wr_out.csv",header=F)
  n_burnin_begin <- round(nrow(wr_out_raw)*pct_burnin_begin)
  n_burnin_end <- round(nrow(wr_out_raw)*pct_burnin_end)
  
  wr_out <- rbind(wr_out,wr_out_raw[n_burnin_begin:n_burnin_end,])
  rm(wr_out_raw)
  
  setwd(results_dir)
}

summary_dir <- paste0("multi_chain_flu_summary_seeds_",min(seeds),"_",max(seeds))

if(!dir.exists(summary_dir)){
  dir.create(summary_dir)
}

setwd(summary_dir)

save(scalars_out,file="scalars_out.RData")
save(wp_out,file="wp_out.RData")
save(tp_out,file="tp_out.RData")
save(dp_out,file="dp_out.RData")
save(wr_out,file="wr_out.RData")
save(model_out,file="model_out.RData")

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

for(i in 2:ncol(scalars_out)){
  print(scalar_plotnames[i-1])
  plot(scalars_out[,i],type="l",main=scalar_plotnames[i-1])
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

print("wp Mean, h3n2")
print(quantile(scalars_out$wp_mean_1,probs=c(0.025,0.5,0.975)))

print("wp Mean, dual")
print(quantile(scalars_out$wp_mean_2,probs=c(0.025,0.5,0.975)))

print("dp Mean, h1n1")
print(quantile(scalars_out$dp_mean_0,probs=c(0.025,0.5,0.975)))

print("dp Mean, h3n2")
print(quantile(scalars_out$dp_mean_1,probs=c(0.025,0.5,0.975)))

print("dp Mean, dual")
print(quantile(scalars_out$dp_mean_2,probs=c(0.025,0.5,0.975)))

print("wr Mean, h1n1")
print(quantile(scalars_out$wr_mean_0,probs=c(0.025,0.5,0.975)))

print("wr Mean, h3n2")
print(quantile(scalars_out$wr_mean_1,probs=c(0.025,0.5,0.975)))

print("wr Mean, dual")
print(quantile(scalars_out$wr_mean_2,probs=c(0.025,0.5,0.975)))

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

print("dp SD, h1n1")
print(quantile(scalars_out$dp_sd_0,probs=c(0.025,0.5,0.975)))

print("wr SD, h1n1")
print(quantile(scalars_out$wr_sd_0,probs=c(0.025,0.5,0.975)))

acp_pr <- function(mcmc_vec){
  return(mean(mcmc_vec[1:(length(mcmc_vec)-1)] != mcmc_vec[2:(length(mcmc_vec))],na.rm=T))
}

scalar_acp_pr <- apply(scalars_out,2,acp_pr)

print("Acceptance Probabilities")
print(scalar_acp_pr)

print("ESS")
print(apply(scalars_out,2,effectiveSize))

sink(file=NULL)


print("Model Traceplots")

pdf(file="Model_traceplots.pdf",width=16,height=16)

par(mfrow=c(5,5))

for(i in 1:ncol(model_out)){
  model_p <- round(table(c(model_out[,i],1,2,3))/(length(model_out[,i])+3),3)

  plot(model_out[,i],main=paste0("Model, ",index_id[i],", p=(",model_p[1],",",model_p[2],",",model_p[3],")"),type="l",xlab="Iteration",ylab="Model")
}

par(mfrow=c(1,1))

dev.off()

## Model estimates
model_est <- apply(model_out,2,function(x){which.max(c(table(c(x,1,2,3))))})
model_est_df <- data.frame(id=individual_data$index,
                           model_est=model_est,
                           subtype=individual_data$subtype)

save(model_est_df,file="model_est.RData")

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
    filter(index==i-1) %>%
    arrange(time)

  model_infer <- which.max(table(c(model_out[,i],1,2,3)))

  which_infer_model <- model_out[,i] == model_infer

  wp_quantile <- quantile(wp_out[which_infer_model,i],probs = c(0.025,0.5,0.975))
  tp_quantile <- quantile(tp_out[which_infer_model,i],probs = c(0.025,0.5,0.975))
  dp_quantile <- quantile(dp_out[which_infer_model,i],probs = c(0.025,0.5,0.975))
  wr_quantile <- quantile(wr_out[which_infer_model,i],probs = c(0.025,0.5,0.975))

  tp_med <- median(tp_out[which_infer_model,i])

  plot_col <- c("blue","red","purple")[individual_data$subtype[i]+1]

  plot(plot_dat$time,plot_dat$viral_load,pch=19,main=paste0("Observed Data, ID ",i-1),col=plot_col,
       xlim=c(min(-wp_quantile[3],plot_dat$time),max(wr_quantile[3],plot_dat$time)),
       ylim=c(0,max(dp_quantile[3],plot_dat$time)))

  if(model_infer != 3){
    lines(c(tp_med-wp_quantile[2],tp_med),c(0,dp_quantile[2]))
    lines(c(tp_med-wp_quantile[1],tp_med),c(0,dp_quantile[1]),lty="dashed")
    lines(c(tp_med-wp_quantile[3],tp_med),c(0,dp_quantile[3]),lty="dashed")
  }

  if(model_infer != 1){
    lines(c(tp_med,tp_med+wr_quantile[2]),c(dp_quantile[2],0))
    lines(c(tp_med,tp_med+wr_quantile[1]),c(dp_quantile[1],0),lty="dashed")
    lines(c(tp_med,tp_med+wr_quantile[3]),c(dp_quantile[3],0),lty="dashed")
  }


}

par(mfrow=c(1,1))
dev.off()

### Summarize the trajectories into subtypes
# t_val <- sort(c(seq(from=-7,to=20,length.out=30),0))
#
# plot_dat_H1N1 <- data.frame(t_val=t_val,
#                             mu_med=rep(NA,length(t_val)),
#                             mu_l=rep(NA,length(t_val)),
#                             mu_h=rep(NA,length(t_val)))
#
#
# for(i in 1:length(t_val)){
#   print(paste0("time: ",t_val[i]))
#   mu_temp <- rep(NA,nrow(scalars_out))
#
#   for(j in 1:nrow(scalars_out)){
#     mu_temp[j] <- ViralLoadRJMCMC::mu(t_val[i],scalars_out$wp_mean_0[j],0,scalars_out$dp_mean_0[j],scalars_out$wr_mean_0[j])
#   }
#
#   plot_dat_H1N1$mu_med[i] <- median(mu_temp)
#   plot_dat_H1N1$mu_l[i] <- quantile(mu_temp,probs=0.025)
#   plot_dat_H1N1$mu_h[i] <- quantile(mu_temp,probs=0.975)
# }
#
# plot_dat_H3N2 <- data.frame(t_val=t_val,
#                             mu_med=rep(NA,length(t_val)),
#                             mu_l=rep(NA,length(t_val)),
#                             mu_h=rep(NA,length(t_val)))
#
#
# for(i in 1:length(t_val)){
#   print(paste0("time: ",t_val[i]))
#   mu_temp <- rep(NA,nrow(scalars_out))
#
#   for(j in 1:nrow(scalars_out)){
#     mu_temp[j] <- ViralLoadRJMCMC::mu(t_val[i],scalars_out$wp_mean_1[j],0,scalars_out$dp_mean_1[j],scalars_out$wr_mean_1[j])
#   }
#
#   plot_dat_H3N2$mu_med[i] <- median(mu_temp)
#   plot_dat_H3N2$mu_l[i] <- quantile(mu_temp,probs=0.025)
#   plot_dat_H3N2$mu_h[i] <- quantile(mu_temp,probs=0.975)
# }
#
# plot_dat_dual <- data.frame(t_val=t_val,
#                             mu_med=rep(NA,length(t_val)),
#                             mu_l=rep(NA,length(t_val)),
#                             mu_h=rep(NA,length(t_val)))
#
#
# for(i in 1:length(t_val)){
#   print(paste0("time: ",t_val[i]))
#   mu_temp <- rep(NA,nrow(scalars_out))
#
#   for(j in 1:nrow(scalars_out)){
#     mu_temp[j] <- ViralLoadRJMCMC::mu(t_val[i],scalars_out$wp_mean_2[j],0,scalars_out$dp_mean_2[j],scalars_out$wr_mean_2[j])
#   }
#
#   plot_dat_dual$mu_med[i] <- median(mu_temp)
#   plot_dat_dual$mu_l[i] <- quantile(mu_temp,probs=0.025)
#   plot_dat_dual$mu_h[i] <- quantile(mu_temp,probs=0.975)
# }
#
# plot_y_min <- lod
# plot_y_max <- min(lod-plot_dat_H1N1$mu_h)
#
# H1N1_trajectory_plot <- plot_dat_H1N1 %>% ggplot(aes(x=t_val,y=lod-mu_med)) +
#   geom_ribbon(aes(ymin=lod-mu_l,ymax=lod-mu_h), alpha=0.5, linewidth = 0, fill="coral") +
#   geom_line(col="coral") +
#   labs(title = "Mean Viral Trajectory, H1N1",fill="Infection Type",color="Infection Type") +
#   xlab("Days since peak") +
#   ylab("Ct") +
#   ylim(plot_y_min,plot_y_max) +
#   theme_linedraw()
#
# print(H1N1_trajectory_plot)

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
end_iter <- (1:(length(seeds)))/length(seeds)*n_iter_total
begin_iter <- c(1,end_iter[1:(length(seeds)-1)]-1)

mcmc_multi_list <- as.mcmc.list(as.mcmc(scalars_out[begin_iter[1]:end_iter[1],c("wp_mean_0","wp_mean_1","wp_mean_2","wp_sd_0","wp_sd_1","wp_sd_2","dp_mean_0","dp_mean_1","dp_mean_2","dp_sd_0","dp_sd_1","dp_sd_2","wr_mean_0","wr_mean_1","wr_mean_2","wr_sd_0","wr_sd_1","wr_sd_2","sigma")]))

for(i in 2:length(seeds)){
  mcmc_multi_list[[i]] <- as.mcmc(scalars_out[begin_iter[i]:end_iter[i],c("wp_mean_0","wp_mean_1","wp_mean_2","wp_sd_0","wp_sd_1","wp_sd_2","dp_mean_0","dp_mean_1","dp_mean_2","dp_sd_0","dp_sd_1","dp_sd_2","wr_mean_0","wr_mean_1","wr_mean_2","wr_sd_0","wr_sd_1","wr_sd_2","sigma")])
}

sink(file="gelman.txt")
print(gelman.diag(mcmc_multi_list,multivariate = F))
sink(file=NULL)


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


png(file="Data_plots_indiv.png",width=1800,height=1800,res=200)
par(mfrow=c(3,3))

id_keep <- c(67,0,24,79,90,94,10,48)

for(i in 1:length(id_keep)){
  plot_dat <- viral_data %>%
    filter(index==id_keep[i]) %>%
    arrange(time_raw)
  
  day_max <- which.max(plot_dat$viral_load)
  
  model_infer <- which.max(table(c(model_out[,(id_keep[i]+1)],1,2,3)))
  
  which_infer_model <- model_out[,(id_keep[i]+1)] == model_infer
  
  wp_quantile <- quantile(wp_out[which_infer_model,(id_keep[i]+1)],probs = c(0.025,0.5,0.975))
  tp_quantile <- quantile(tp_out[which_infer_model,(id_keep[i]+1)],probs = c(0.025,0.5,0.975))
  dp_quantile <- quantile(dp_out[which_infer_model,(id_keep[i]+1)],probs = c(0.025,0.5,0.975))
  wr_quantile <- quantile(wr_out[which_infer_model,(id_keep[i]+1)],probs = c(0.025,0.5,0.975))
  
  tp_med <- median(tp_out[which_infer_model,(id_keep[i]+1)]) + day_max
  
  plot_col <- c("blue","red","purple")[individual_data$subtype[(id_keep[i]+1)]+1]
  plot_subtype <- c("H1N1","H3N2","Dual")[individual_data$subtype[(id_keep[i]+1)]+1]
  
  if(model_infer == 1){
    limx = c(min(tp_med-wp_quantile[3],0),max(tp_med,plot_dat$time_raw))
  } else if(model_infer == 2){
    limx = c(min(tp_med-wp_quantile[3],0),max(tp_med+wr_quantile[3],plot_dat$time_raw))
  } else if(model_infer == 3){
    limx = c(min(tp_med,0),max(tp_med+wr_quantile[3],plot_dat$time_raw))
  } else{
    UNTITLED()
  }
  
  plot(plot_dat$time_raw,settings$lod-plot_dat$viral_load,
       pch=ifelse(plot_dat$viral_load > 0,19,1),
       main=paste0("Observed Data, ID ",id_keep[i]," (",plot_subtype,")"),
       col=plot_col,
       xlim=limx,
       ylim=c(settings$lod,18),
       xlab="Day of Fair",
       ylab="Ct")
  
  if(model_infer != 3){
    lines(c(tp_med-wp_quantile[2],tp_med),c(settings$lod,settings$lod-dp_quantile[2]))
    lines(c(tp_med-wp_quantile[1],tp_med),c(settings$lod,settings$lod-dp_quantile[1]),lty="dashed")
    lines(c(tp_med-wp_quantile[3],tp_med),c(settings$lod,settings$lod-dp_quantile[3]),lty="dashed")
  }
  
  if(model_infer != 1){
    lines(c(tp_med,tp_med+wr_quantile[2]),c(settings$lod-dp_quantile[2],settings$lod))
    lines(c(tp_med,tp_med+wr_quantile[1]),c(settings$lod-dp_quantile[1],settings$lod),lty="dashed")
    lines(c(tp_med,tp_med+wr_quantile[3]),c(settings$lod-dp_quantile[3],settings$lod),lty="dashed")
  }
  
  abline(v=max(plot_dat$time_raw),lty="dashed",col="grey50")
  
}

par(mfrow=c(1,1))
dev.off()

png(file="data_example.png",width=1800,height=600,res=200)
par(mfrow=c(1,3))

which_example <- c(2,0,99)
example_subtype <- c("H1N1","H3N2","Dual")

for(i in 1:length(which_example)){
  plot_dat <- viral_data %>%
    filter(index==which_example[i]) %>%
    arrange(time_raw)
  
  day_max <- which.max(plot_dat$viral_load)
  
  model_infer <- which.max(table(c(model_out[,(which_example[i]+1)],1,2,3)))
  
  which_infer_model <- model_out[,(which_example[i]+1)] == model_infer
  
  wp_quantile <- quantile(wp_out[which_infer_model,(which_example[i]+1)],probs = c(0.025,0.5,0.975))
  tp_quantile <- quantile(tp_out[which_infer_model,(which_example[i]+1)],probs = c(0.025,0.5,0.975))
  dp_quantile <- quantile(dp_out[which_infer_model,(which_example[i]+1)],probs = c(0.025,0.5,0.975))
  wr_quantile <- quantile(wr_out[which_infer_model,(which_example[i]+1)],probs = c(0.025,0.5,0.975))
  
  tp_med <- median(tp_out[which_infer_model,(which_example[i]+1)]) + day_max
  
  plot_col <- c("blue","red","purple")[individual_data$subtype[(which_example[i]+1)]+1]
  
  # if(model_infer == 1){
  #   limx = c(min(tp_med-wp_quantile[3],0),max(tp_med,plot_dat$time_raw))
  # } else if(model_infer == 2){
  #   limx = c(min(tp_med-wp_quantile[3],0),max(tp_med+wr_quantile[3],plot_dat$time_raw))
  # } else if(model_infer == 3){
  #   limx = c(min(tp_med,0),max(tp_med+wr_quantile[3],plot_dat$time_raw))
  # } else{
  #   UNTITLED()
  # }
  
  limx <- c(0,7)
  
  plot(plot_dat$time_raw,settings$lod-plot_dat$viral_load,pch=ifelse(plot_dat$viral_load==0,1,19),main=paste0("Observed Data, ID ",which_example[i]," (",example_subtype[i],")"),col=plot_col,
       cex=1.5,
       xlim=limx,
       ylim=c(settings$lod,18),
       xlab="Day of Fair",
       ylab="Ct")
  
  # if(model_infer != 3){
  #   lines(c(tp_med-wp_quantile[2],tp_med),c(settings$lod,settings$lod-dp_quantile[2]))
  #   lines(c(tp_med-wp_quantile[1],tp_med),c(settings$lod,settings$lod-dp_quantile[1]),lty="dashed")
  #   lines(c(tp_med-wp_quantile[3],tp_med),c(settings$lod,settings$lod-dp_quantile[3]),lty="dashed")
  # }
  # 
  # if(model_infer != 1){
  #   lines(c(tp_med,tp_med+wr_quantile[2]),c(settings$lod-dp_quantile[2],settings$lod))
  #   lines(c(tp_med,tp_med+wr_quantile[1]),c(settings$lod-dp_quantile[1],settings$lod),lty="dashed")
  #   lines(c(tp_med,tp_med+wr_quantile[3]),c(settings$lod-dp_quantile[3],settings$lod),lty="dashed")
  # }
  
  # abline(v=max(plot_dat$time_raw),lty="dashed",col="grey50")
  
}

par(mfrow=c(1,1))
dev.off()

### Violin plots
### Parameter violin plots
wp_mean_violin_data <- data.frame(subtype=factor(c(rep("H1N1",nrow(scalars_out)),rep("H3N2",nrow(scalars_out)),rep("Dual",nrow(scalars_out))),levels=c("H1N1","H3N2","Dual")),
                                  value=c(scalars_out$wp_mean_0,scalars_out$wp_mean_1,scalars_out$wp_mean_2))

wp_mean_violin <- ggplot(data=wp_mean_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=2.5) +
  scale_fill_manual(values=c("skyblue2","coral","purple3"),guide="none") +
  labs(x="",y="Mean (Days)",fill="Subtype") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
  ggtitle("Mean Proliferation Duration")

wp_mean_violin

dp_mean_violin_data <- data.frame(subtype=factor(c(rep("H1N1",nrow(scalars_out)),rep("H3N2",nrow(scalars_out)),rep("Dual",nrow(scalars_out))),levels=c("H1N1","H3N2","Dual")),
                                  value=c(scalars_out$dp_mean_0,scalars_out$dp_mean_1,scalars_out$dp_mean_2))

dp_mean_violin <- ggplot(data=dp_mean_violin_data,aes(x=factor(subtype),y=45-value,fill=subtype,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=4) +
  scale_fill_manual(values=c("skyblue2","coral","purple3"),guide="none") +
  scale_y_reverse() +
  labs(x="",y="Mean Peak Ct",fill="Subtype") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
  ggtitle("Mean Viral Load Peak")

dp_mean_violin

wr_mean_violin_data <- data.frame(subtype=factor(c(rep("H1N1",nrow(scalars_out)),rep("H3N2",nrow(scalars_out)),rep("Dual",nrow(scalars_out))),levels=c("H1N1","H3N2","Dual")),
                                  value=c(scalars_out$wr_mean_0,scalars_out$wr_mean_1,scalars_out$wr_mean_2))

wr_mean_violin <- ggplot(data=wr_mean_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=4) +
  scale_fill_manual(values=c("skyblue2","coral","purple3"),guide="none") +
  labs(x="",y="Mean (Days)",fill="Subtype") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
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
    ggtitle("SD, Proliferation Duration")
} else{
  wp_sd_violin_data <- data.frame(subtype=factor(c(rep("H1N1",nrow(scalars_out)),rep("H3N2",nrow(scalars_out)),rep("Dual",nrow(scalars_out))),levels=c("H1N1","H3N2","Dual")),
                                  value=c(scalars_out$wp_sd_0,scalars_out$wp_sd_1,scalars_out$wp_sd_2))
  
  wp_sd_violin <- ggplot(data=wp_sd_violin_data,aes(x=factor(subtype),y=value,fill=subtype,alpha = 0.8)) +
    geom_violin(trim=F,linewidth=0.2,adjust=2.5) +
    scale_fill_manual(values=c("skyblue2","coral","purple3"),guide="none") +
    labs(x="",y="SD Proliferation",fill="Subtype") +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size=11)) +
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
    ggtitle("SD, Viral Load Peak")
} else{
  dp_sd_violin_data <- data.frame(subtype=factor(c(rep("H1N1",nrow(scalars_out)),rep("H3N2",nrow(scalars_out)),rep("Dual",nrow(scalars_out))),levels=c("H1N1","H3N2","Dual")),
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
    ggtitle("SD, Clearance Duration")
} else{
  wr_sd_violin_data <- data.frame(subtype=factor(c(rep("H1N1",nrow(scalars_out)),rep("H3N2",nrow(scalars_out)),rep("Dual",nrow(scalars_out))),levels=c("H1N1","H3N2","Dual")),
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

posterior_violin_plots_trajectories <- ggarrange(wp_mean_violin,
                                                 dp_mean_violin,
                                                 wr_mean_violin,
                                                 wp_sd_violin,
                                                 dp_sd_violin,
                                                 wr_sd_violin,
                                                 nrow=2,ncol=3)

annotate_figure(posterior_violin_plots_trajectories, top = text_grob("Posterior Distributions", 
                                                                     face = "bold", size = 14))

png(filename="posterior_violins_mean_sd_multichain.png",width=1400,height=1000,res=200)
print(posterior_violin_plots_trajectories)
dev.off()


### Parameter difference violin plots
wp_mean_diff_violin_data <- data.frame(subtype=factor(c(rep("H1N1-H3N2",nrow(scalars_out)),rep("Dual-H1N1",nrow(scalars_out)),rep("Dual-H3N2",nrow(scalars_out))),levels=c("H1N1-H3N2","Dual-H1N1","Dual-H3N2")),
                                       value=c(scalars_out$wp_mean_0-scalars_out$wp_mean_1,scalars_out$wp_mean_2 - scalars_out$wp_mean_0,scalars_out$wp_mean_2 - scalars_out$wp_mean_1))

wp_mean_diff_violin <- ggplot(data=wp_mean_diff_violin_data,aes(x=factor(subtype),y=value,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=2.5,fill="grey30",width=0.65) +
  geom_hline(yintercept=0) +
  # scale_fill_manual(values=c("coral","skyblue2","purple3"),guide="none") +
  labs(x="",y="Mean Proliferation Difference") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle("Mean Proliferation, Difference")

wp_mean_diff_violin

dp_mean_diff_violin_data <- data.frame(subtype=factor(c(rep("H1N1-H3N2",nrow(scalars_out)),rep("Dual-H1N1",nrow(scalars_out)),rep("Dual-H3N2",nrow(scalars_out))),levels=c("H1N1-H3N2","Dual-H1N1","Dual-H3N2")),
                                       value=c(scalars_out$dp_mean_0-scalars_out$dp_mean_1,scalars_out$dp_mean_2-scalars_out$dp_mean_0,scalars_out$dp_mean_2-scalars_out$dp_mean_1))

dp_mean_diff_violin <- ggplot(data=dp_mean_diff_violin_data,aes(x=factor(subtype),y=-value,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=2.5,fill="grey30",width=0.65) +
  geom_hline(yintercept=0) +
  # scale_fill_manual(values=c("coral","skyblue2","purple3"),guide="none") +
  labs(x="",y="Mean Peak Ct Difference") +
  scale_y_reverse() +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle("Mean Peak Ct, Difference")

dp_mean_diff_violin

wr_mean_diff_violin_data <- data.frame(subtype=factor(c(rep("H1N1-H3N2",nrow(scalars_out)),rep("Dual-H1N1",nrow(scalars_out)),rep("Dual-H3N2",nrow(scalars_out))),levels=c("H1N1-H3N2","Dual-H1N1","Dual-H3N2")),
                                       value=c(scalars_out$wr_mean_0-scalars_out$wr_mean_1,scalars_out$wr_mean_2-scalars_out$wr_mean_0,scalars_out$wr_mean_2-scalars_out$wr_mean_1))

wr_mean_diff_violin <- ggplot(data=wr_mean_diff_violin_data,aes(x=factor(subtype),y=value,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=4,fill="grey30",width=0.65) +
  geom_hline(yintercept=0) +
  # scale_fill_manual(values=c("coral","skyblue2","purple3"),guide="none") +
  labs(x="",y="Mean Clearance Difference") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle("Mean Clearance, Difference")

wr_mean_diff_violin

wp_sd_diff_violin_data <- data.frame(subtype=factor(c(rep("H1N1-H3N2",nrow(scalars_out)),rep("Dual-H1N1",nrow(scalars_out)),rep("Dual-H3N2",nrow(scalars_out))),levels=c("H1N1-H3N2","Dual-H1N1","Dual-H3N2")),
                                     value=c(scalars_out$wp_sd_0-scalars_out$wp_sd_1,scalars_out$wp_sd_2-scalars_out$wp_sd_0,scalars_out$wp_sd_2-scalars_out$wp_sd_1))

wp_sd_diff_violin <- ggplot(data=wp_sd_diff_violin_data,aes(x=factor(subtype),y=value,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=2.5,fill="grey30") +
  # scale_fill_manual(values=c("coral","skyblue2","purple3"),guide="none") +
  labs(x="",y="SD proliferation") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle("SD, Proliferation, Difference")

wp_sd_diff_violin

dp_sd_diff_violin_data <- data.frame(subtype=factor(c(rep("H1N1-H3N2",nrow(scalars_out)),rep("Dual-H1N1",nrow(scalars_out)),rep("Dual-H3N2",nrow(scalars_out))),levels=c("H1N1-H3N2","Dual-H1N1","Dual-H3N2")),
                                     value=c(scalars_out$dp_sd_0-scalars_out$dp_sd_1,scalars_out$dp_sd_2-scalars_out$dp_sd_0,scalars_out$dp_sd_2-scalars_out$dp_sd_1))

dp_sd_diff_violin <- ggplot(data=dp_sd_diff_violin_data,aes(x=factor(subtype),y=value,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=2.5,fill="grey30") +
  # scale_fill_manual(values=c("coral","skyblue2","purple3"),guide="none") +
  labs(x="",y=expression(d[p])) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle("SD, Viral Load Peak,\n
          Difference")

dp_sd_diff_violin

wr_sd_diff_violin_data <- data.frame(subtype=factor(c(rep("H1N1-H3N2",nrow(scalars_out)),rep("Dual-H1N1",nrow(scalars_out)),rep("Dual-H3N2",nrow(scalars_out))),levels=c("H1N1-H3N2","Dual-H1N1","Dual-H3N2")),
                                     value=c(scalars_out$wr_sd_0-scalars_out$wr_sd_1,scalars_out$wr_sd_2-scalars_out$wr_sd_0,scalars_out$wr_sd_2-scalars_out$wr_sd_1))

wr_sd_diff_violin <- ggplot(data=wr_sd_diff_violin_data,aes(x=factor(subtype),y=value,alpha = 0.8)) +
  geom_violin(trim=F,linewidth=0.2,adjust=2.5,fill="grey30") +
  # scale_fill_manual(values=c("coral","skyblue2","purple3"),guide="none") +
  labs(x="",y=expression(omega[r])) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ggtitle("SD, Clearance Duration,\n
          Difference")

wr_sd_diff_violin

posterior_diff_violin_plots_trajectories <- ggarrange(wp_mean_diff_violin,
                                                      dp_mean_diff_violin,
                                                      wr_mean_diff_violin,
                                                      # wp_sd_diff_violin,
                                                      # dp_sd_diff_violin,
                                                      # wr_sd_diff_violin,
                                                      nrow=3,ncol=1)

posterior_diff_violin_plots_trajectories_horizontal <- ggarrange(wp_mean_diff_violin,
                                                                 dp_mean_diff_violin,
                                                                 wr_mean_diff_violin,
                                                                 # wp_sd_diff_violin,
                                                                 # dp_sd_diff_violin,
                                                                 # wr_sd_diff_violin,
                                                                 nrow=1,ncol=3)

annotate_figure(posterior_diff_violin_plots_trajectories, top = text_grob("Posterior Differences", 
                                                                          face = "bold", size = 14))

annotate_figure(posterior_diff_violin_plots_trajectories_horizontal, top = text_grob("Posterior Differences", 
                                                                                     face = "bold", size = 14))

png(filename="posterior_diff_violins_mean_sd_multichain.png",width=1000,height=2000,res=200)
print(posterior_diff_violin_plots_trajectories)
dev.off()

png(filename="posterior_diff_violins_mean_sd_multichain_horizontal.png",width=2100,height=900,res=200)
print(posterior_diff_violin_plots_trajectories_horizontal)
dev.off()