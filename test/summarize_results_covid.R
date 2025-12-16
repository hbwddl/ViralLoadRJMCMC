library(Rcpp)
library(dplyr)
library(truncnorm)
library(ggpubr)
library(ggplot2)

dec.precision <- 6

pct_burnin_begin <- 0.5
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

scalar_plotnames <- c("wp Mean, Asymptomatic",
                      "wp Mean, Symptomatic",
                      "wp SD, Asymptomatic",
                      "wp SD, Symptomatic",
                      "Mean Peak Ct, Asymptomatic",
                      "Mean Peak Ct, Symptomatic",
                      "dp SD, Asymptomatic",
                      "dp SD, Symptomatic",
                      "tp SD, Asymptomatic",
                      "tp SD, Symptomatic",
                      "wr Mean, Asymptomatic",
                      "wr Mean, Symptomatic",
                      "wr SD, Asymptomatic",
                      "wr SD, Symptomatic",
                      "Sigma",
                      "Log Likelihood")

print("Scalar plots")

pdf(file="Scalar_Plots.pdf",width=8,height=8)
par(mfrow=c(2,2))

for(i in 2:ncol(scalars_out)){
  print(scalar_plotnames[i-1])
  if(i == 6 | i==7){
    plot(settings$lod-scalars_out[,i],type="l",main=scalar_plotnames[i-1])
  } else{
    plot(scalars_out[,i],type="l",main=scalar_plotnames[i-1])
  }
}

par(mfrow=c(1,1))
dev.off()

print("Scalar Difference Plots")

pdf(file="Scalar_Diff_Plots.pdf",width=14,height=5)

par(mfrow=c(1,3))
plot(scalars_out$wp_mean_1 - scalars_out$wp_mean_0,type="l",main="WP Mean Difference, Symptomatic-Asymptomatic",xlab="Iteration",ylab="Difference")
abline(h=0)

plot(scalars_out$dp_mean_1 - scalars_out$dp_mean_0,type="l",main="DP Mean Difference, Symptomatic-Asymptomatic",xlab="Iteration",ylab="Difference")
abline(h=0)

plot(scalars_out$wr_mean_1 - scalars_out$wr_mean_0,type="l",main="WR Mean Difference, Symptomatic-Asymptomatic",xlab="Iteration",ylab="Difference")
abline(h=0)

par(mfrow=c(1,1))

dev.off()

print("Quantiles")

sink(file="Quantiles.txt")
print("wp Mean, Asymptomatic")
print(quantile(scalars_out$wp_mean_0,probs=c(0.025,0.5,0.975)))

print("wp Mean, Symptomatic")
print(quantile(scalars_out$wp_mean_1,probs=c(0.025,0.5,0.975)))

print("Peak Ct, Asymptomatic")
print(quantile(settings$lod - scalars_out$dp_mean_0,probs=c(0.025,0.5,0.975)))

print("Peak Ct, Symptomatic")
print(quantile(settings$lod - scalars_out$dp_mean_1,probs=c(0.025,0.5,0.975)))

print("wr Mean, Asymptomatic")
print(quantile(scalars_out$wr_mean_0,probs=c(0.025,0.5,0.975)))

print("wr Mean, Symptomatic")
print(quantile(scalars_out$wr_mean_1,probs=c(0.025,0.5,0.975)))

print("WP Mean Difference, Symptomatic-Asymptomatic")
print(quantile(scalars_out$wp_mean_1 - scalars_out$wp_mean_0,probs=c(0.025,0.5,0.975)))

print("DP Mean Difference, Symptomatic-Asymptomatic")
print(quantile(scalars_out$dp_mean_1 - scalars_out$dp_mean_0,probs=c(0.025,0.5,0.975)))

print("WR Mean Difference, Symptomatic-Asymptomatic")
print(quantile(scalars_out$wr_mean_1 - scalars_out$wr_mean_0,probs=c(0.025,0.5,0.975)))


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

## Forest plot
forest_plot_dat <- data.frame(analysis=rep(c("Kissler","Kissler","RJMCMC","RJMCMC"),3),
                              parameter=rep(c("Proliferation","Peak Ct","Clearance"),each=4),
                              symptomatic=rep(c("Symptomatic","Asymptomatic"),6),
                              median=c(3.3,3.4,
                                       quantile(scalars_out$wp_mean_1,probs=0.5),quantile(scalars_out$wp_mean_0,probs=0.5),
                                       22.2,22.4,
                                       quantile(settings$lod - scalars_out$dp_mean_1,probs=0.5),quantile(settings$lod - scalars_out$dp_mean_0,probs=0.5),
                                       10.9,7.8,
                                       quantile(scalars_out$wr_mean_1,probs=0.5),quantile(scalars_out$wr_mean_0,probs=0.5)),
                              lower=c(1.9,2.5,
                                      quantile(scalars_out$wp_mean_1,probs=0.025),quantile(scalars_out$wp_mean_0,probs=0.025),
                                      19.1,20.2,
                                      quantile(settings$lod - scalars_out$dp_mean_1,probs=0.025),quantile(settings$lod - scalars_out$dp_mean_0,probs=0.025),
                                      7.8,6.1,
                                      quantile(scalars_out$wr_mean_1,probs=0.025),quantile(scalars_out$wr_mean_0,probs=0.025)),
                              upper=c(5.1,4.5,
                                      quantile(scalars_out$wp_mean_1,probs=0.975),quantile(scalars_out$wp_mean_0,probs=0.975),
                                      25.0,24.5,
                                      quantile(settings$lod - scalars_out$dp_mean_1,probs=0.975),quantile(settings$lod - scalars_out$dp_mean_0,probs=0.975),
                                      14.2,9.7,
                                      quantile(scalars_out$wr_mean_1,probs=0.975),quantile(scalars_out$wr_mean_0,probs=0.975)))

t_min <- -8
t_max <- 15

plot_dat_symptomatic <- data.frame(t_val = c(t_min,
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
plot_y_max <- min(settings$lod-max(plot_dat_symptomatic$h_val),19.1)

symptomatic_trajectory_plot <- plot_dat_symptomatic %>% ggplot(aes(x=t_val,y=settings$lod-med_val)) +
  geom_ribbon(aes(ymin=settings$lod-l_val,ymax=settings$lod-h_val), alpha=0.6, linewidth = 0, fill="coral") +
  geom_line(col="coral") +
  labs(title = "Mean Viral Trajectory, Symptomatic",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Ct") +
  ylim(plot_y_min,plot_y_max) +
  annotate("segment",x=-5.10000,y=settings$lod,xend=0,yend=19.1,linetype="dotted",alpha=0.6) +
  annotate("segment",x=-3.3,y=settings$lod,xend=0,yend=22.2,linetype="dashed",alpha=0.6) +
  annotate("segment",x=-1.9,y=settings$lod,xend=0,yend=25,linetype="dotted",alpha=0.6) +
  annotate("segment",x=0,y=19.1,xend=14.2,yend=settings$lod,linetype="dotted",alpha=0.6) +
  annotate("segment",x=0,y=22.2,xend=10.9,yend=settings$lod,linetype="dashed",alpha=0.6) +
  annotate("segment",x=0,y=25,xend=7.8,yend=settings$lod,linetype="dotted",alpha=0.6) +
  theme_linedraw()

print(symptomatic_trajectory_plot)                        


plot_dat_asymptomatic <- data.frame(t_val = c(t_min,
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

plot_y_min <- settings$lod
plot_y_max <- min(settings$lod-max(plot_dat_asymptomatic$h_val),19.1)

asymptomatic_trajectory_plot <- plot_dat_asymptomatic %>% ggplot(aes(x=t_val,y=settings$lod-med_val)) +
  geom_ribbon(aes(ymin=settings$lod-l_val,ymax=settings$lod-h_val), alpha=0.6, linewidth = 0, fill="skyblue4") +
  geom_line(col="skyblue3") +
  labs(title = "Mean Viral Trajectory, Asymptomatic",fill="Infection Type",color="Infection Type") +
  xlab("Days since peak") +
  ylab("Ct") +
  ylim(plot_y_min,plot_y_max) +
  annotate("segment",x=-4.50000,y=settings$lod,xend=0,yend=20.2,linetype="dotted",alpha=0.6) +
  annotate("segment",x=-3.4,y=settings$lod,xend=0,yend=22.4,linetype="dashed",alpha=0.6) +
  annotate("segment",x=-2.5,y=settings$lod,xend=0,yend=24.5,linetype="dotted",alpha=0.6) +
  annotate("segment",x=0,y=20.2,xend=9.7,yend=settings$lod,linetype="dotted",alpha=0.6) +
  annotate("segment",x=0,y=22.4,xend=7.8,yend=settings$lod,linetype="dashed",alpha=0.6) +
  annotate("segment",x=0,y=24.5,xend=6.1,yend=settings$lod,linetype="dotted",alpha=0.6) +
  theme_linedraw()

print(asymptomatic_trajectory_plot)   

print("Trajectory Plots")

png(filename="trajectory_plots.png",width=2000,height=700,res=200)

print(ggarrange(symptomatic_trajectory_plot,asymptomatic_trajectory_plot,ncol=2))

dev.off()

## Estimated models
## Model estimates
model_est <- apply(model_out,2,function(x){which.max(c(table(c(x,1,2,3))))})
model_est_df <- data.frame(id=individual_data$index,
                           model_est=model_est,
                           subtype=individual_data$subtype)

save(model_est_df,file="model_est.RData")