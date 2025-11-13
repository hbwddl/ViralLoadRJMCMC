library(Rcpp)
library(dplyr)
library(truncnorm)
library(ggpubr)
library(ggplot2)
library(coda)

dec.precision <- 6

pct_burnin_begin <- 0.6
pct_burnin_end <- 0.99

results_dir <- "~/Documents/Research/Within-Host/RJMCMC_Results"

setwd(results_dir)

seeds <- c(1111:1115)

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
# lod <- settings$lod
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
  geom_ribbon(aes(ymin=lod-l_val,ymax=lod-h_val), alpha=0.5, linewidth = 0, fill="coral") +
  geom_line(col="coral") +
  labs(title = "Mean Viral Trajectory, H1N1",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Ct") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

print(H1N1_trajectory_plot)

H3N2_trajectory_plot <- plot_dat_H3N2 %>% ggplot(aes(x=t_val,y=lod-med_val)) +
  geom_ribbon(aes(ymin=lod-l_val,ymax=lod-h_val), alpha=0.5, linewidth = 0, fill="skyblue2") +
  geom_line(col="skyblue2") +
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