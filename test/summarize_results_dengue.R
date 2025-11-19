library(Rcpp)
library(dplyr)
library(truncnorm)
library(ggpubr)
library(ggplot2)
library(coda)

dec.precision <- 6

pct_burnin_begin <- 0.6
pct_burnin_end <- 0.99

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

scalar_plotnames <- c("wp Mean, No leakage",
                      "wp Mean, Leakage",
                      "wp SD, No leakage",
                      "wp SD, Leakage",
                      "dp Mean, No leakage",
                      "dp Mean, Leakage",
                      "dp SD, No leakage",
                      "dp SD, Leakage",
                      "tp SD, No leakage",
                      "tp SD, Leakage",
                      "wr Mean, No leakage",
                      "wr Mean, Leakage",
                      "wr SD, No leakage",
                      "wr SD, Leakage",
                      "Sigma",
                      "Log Likelihood")

print("Scalar plots")

pdf(file="Scalar_Plots.pdf",width=8,height=8)
par(mfrow=c(2,2))

for(i in 2:ncol(scalars_out)){
  print(scalar_plotnames[i-1])
  plot(scalars_out[,i],type="l",main=scalar_plotnames[i-1])
}

par(mfrow=c(1,1))
dev.off()

print("Scalar Difference Plots")

pdf(file="Scalar_Diff_Plots.pdf",width=14,height=5)

par(mfrow=c(1,3))
plot(scalars_out$wp_mean_1 - scalars_out$wp_mean_0,type="l",main="WP Mean Difference, Leakage-No Leakage",xlab="Iteration",ylab="Difference")
abline(h=0)

plot(scalars_out$dp_mean_1 - scalars_out$dp_mean_0,type="l",main="DP Mean Difference, Leakage-No Leakage",xlab="Iteration",ylab="Difference")
abline(h=0)

plot(scalars_out$wr_mean_1 - scalars_out$wr_mean_0,type="l",main="WR Mean Difference, Leakage-No Leakage",xlab="Iteration",ylab="Difference")
abline(h=0)

par(mfrow=c(1,1))

dev.off()

print("Quantiles")

sink(file="Quantiles.txt")
print("wp Mean, No Leakage")
print(quantile(scalars_out$wp_mean_0,probs=c(0.025,0.5,0.975)))

print("wp Mean, Leakage")
print(quantile(scalars_out$wp_mean_1,probs=c(0.025,0.5,0.975)))

print("dp Mean, No Leakage")
print(quantile(scalars_out$dp_mean_0,probs=c(0.025,0.5,0.975)))

print("dp Mean, Leakage")
print(quantile(scalars_out$dp_mean_1,probs=c(0.025,0.5,0.975)))

print("wr Mean, No Leakage")
print(quantile(scalars_out$wr_mean_0,probs=c(0.025,0.5,0.975)))

print("wr Mean, Leakage")
print(quantile(scalars_out$wr_mean_1,probs=c(0.025,0.5,0.975)))

print("WP Mean Difference, Leakage-No Leakage")
print(quantile(scalars_out$wp_mean_1 - scalars_out$wp_mean_0,probs=c(0.025,0.5,0.975)))

print("DP Mean Difference, Leakage-No Leakage")
print(quantile(scalars_out$dp_mean_1 - scalars_out$dp_mean_0,probs=c(0.025,0.5,0.975)))

print("WR Mean Difference, Leakage-No Leakage")
print(quantile(scalars_out$wr_mean_1 - scalars_out$wr_mean_0,probs=c(0.025,0.5,0.975)))


acp_pr <- function(mcmc_vec){
  return(mean(mcmc_vec[1:(length(mcmc_vec)-1)] != mcmc_vec[2:(length(mcmc_vec))],na.rm=T))
}

scalar_acp_pr <- apply(scalars_out,2,acp_pr)

print("Acceptance Probabilities")
print(scalar_acp_pr)

print("ESS")
print(apply(scalars_out,2,coda::effectiveSize))

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
    filter(index==i-1) %>%
    arrange(time)

  wp_quantile <- quantile(wp_out[,i],probs = c(0.025,0.5,0.975))
  tp_quantile <- quantile(tp_out[,i],probs = c(0.025,0.5,0.975))
  dp_quantile <- quantile(dp_out[,i],probs = c(0.025,0.5,0.975))
  wr_quantile <- quantile(wr_out[,i],probs = c(0.025,0.5,0.975))

  tp_med <- median(tp_out[,i])

  model_infer <- which.max(table(c(model_out[,i],1,2,3)))

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

t_min <- -20
t_max <- 15

plot_dat_leakage <- data.frame(t_val = c(t_min,
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

plot_y_min <- settings$lod
plot_y_max <- min(settings$lod-max(plot_dat_leakage$h_val,plot_dat_no_leakage$h_val))

leakage_trajectory_plot <- plot_dat_leakage %>% ggplot(aes(x=t_val,y=settings$lod-med_val)) +
  geom_ribbon(aes(ymin=settings$lod-l_val,ymax=settings$lod-h_val), alpha=0.6, linewidth = 0, fill="coral") +
  geom_line(col="coral") +
  labs(title = "Mean Viral Trajectory, Plasma Leakage",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Viral Count (Log 10)") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

print(leakage_trajectory_plot)                        


plot_dat_no_leakage <- data.frame(t_val = c(t_min,
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

no_leakage_trajectory_plot <- plot_dat_no_leakage %>% ggplot(aes(x=t_val,y=settings$lod-med_val)) +
  geom_ribbon(aes(ymin=settings$lod-l_val,ymax=settings$lod-h_val), alpha=0.6, linewidth = 0, fill="skyblue4") +
  geom_line(col="skyblue3") +
  labs(title = "Mean Viral Trajectory, No Leakage",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Viral Count (Log 10)") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

print(no_leakage_trajectory_plot)   

print("Trajectory Plots Full")

png(filename="trajectory_plots_full.png",width=2000,height=700,res=200)

print(ggarrange(leakage_trajectory_plot,no_leakage_trajectory_plot,ncol=2))

dev.off()

t_min <- -1
t_max <- 5

plot_y_min <- settings$lod
plot_y_max <- min(settings$lod-max(plot_dat_leakage$h_val,plot_dat_no_leakage$h_val))

half_dat_leakage <- data.frame(t_val = c(0,
                                         quantile(scalars_out$wr_mean_1,probs=0.025),
                                         quantile(scalars_out$wr_mean_1,probs=0.5),
                                         quantile(scalars_out$wr_mean_1,probs=0.975),
                                         t_max),
                               med_val = c(quantile(scalars_out$dp_mean_1,probs=0.5),
                                           mu(quantile(scalars_out$wr_mean_1,probs=0.025),
                                              quantile(scalars_out$wp_mean_1,probs=0.5),
                                              0,
                                              quantile(scalars_out$dp_mean_1,probs=0.5),
                                              quantile(scalars_out$wr_mean_1,probs=0.5)),
                                           0,
                                           0,
                                           0),
                               l_val = c(quantile(scalars_out$dp_mean_1,probs=0.025),
                                         0,
                                         0,
                                         0,
                                         0),
                               h_val = c(quantile(scalars_out$dp_mean_1,probs=0.975),
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

plot_y_min <- 0
plot_y_max <- settings$lod

half_leakage_trajectory_plot <- half_dat_leakage %>% ggplot(aes(x=t_val,y=med_val)) +
  geom_ribbon(aes(ymin=l_val,ymax=h_val), alpha=0.6, linewidth = 0, fill="coral") +
  geom_line(col="coral") +
  labs(title = "Mean Clearance Trajectory, Plasma Leakage",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Viral Count (Log 10)") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

print(half_leakage_trajectory_plot)                        

half_dat_no_leakage <- data.frame(t_val = c(0,
                                            quantile(scalars_out$wr_mean_0,probs=0.025),
                                            quantile(scalars_out$wr_mean_0,probs=0.5),
                                            quantile(scalars_out$wr_mean_0,probs=0.975),
                                            t_max),
                                  med_val = c(quantile(scalars_out$dp_mean_0,probs=0.5),
                                              mu(quantile(scalars_out$wr_mean_0,probs=0.025),
                                                 quantile(scalars_out$wp_mean_0,probs=0.5),
                                                 0,
                                                 quantile(scalars_out$dp_mean_0,probs=0.5),
                                                 quantile(scalars_out$wr_mean_0,probs=0.5)),
                                              0,
                                              0,
                                              0),
                                  l_val = c(quantile(scalars_out$dp_mean_0,probs=0.025),
                                            0,
                                            0,
                                            0,
                                            0),
                                  h_val = c(quantile(scalars_out$dp_mean_0,probs=0.975),
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

half_no_leakage_trajectory_plot <- half_dat_no_leakage %>% ggplot(aes(x=t_val,y=med_val)) +
  geom_ribbon(aes(ymin=l_val,ymax=h_val), alpha=0.6, linewidth = 0, fill="skyblue4") +
  geom_line(col="skyblue3") +
  labs(title = "Mean Clearance Trajectory, No Leakage",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Viral Count (Log 10)") +
  ylim(plot_y_min,plot_y_max) +
  theme_linedraw()

print(half_no_leakage_trajectory_plot)   

print("Trajectory Plots")

png(filename="trajectory_plots.png",width=1000,height=1200,res=200)

print(ggarrange(half_leakage_trajectory_plot,half_no_leakage_trajectory_plot,ncol=1))

dev.off()

