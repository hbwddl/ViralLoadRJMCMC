### Log likelihood for toy data
library(ViralLoadRJMCMC)
set.seed(101)

sink(file="lh_test_r.txt")
wp_data <- c(3,4,2,3,4,5)
dp_data <- c(35,34,33,20,40,30)
tp_data <- c(1,0.5,0.1,-0.5,-1.5,-1)
wr_data <- c(5,4,6,7,3,4)

subtype <- c(0,1,0,1,0,1)

sigma <- 5

model <- c(1,1,2,2,3,3)

sensitivity <- 0.99

obs_data <- data.frame(index=c(rep(0,5),
                               rep(1,4),
                               rep(2,6),
                               rep(3,6),
                               rep(4,6),
                               rep(5,5)),
                       time=c(-4:0,
                              -3:0,
                              -2:3,
                              -3:2,
                              -1:4,
                              0:4))

obs_data$mu <- NA

for(i in 1:nrow(obs_data)){
  obs_data$mu[i] <- mu(obs_data$time[i],
                               wp_data[obs_data$index[i]+1],
                               tp_data[obs_data$index[i]+1],
                               dp_data[obs_data$index[i]+1],
                               wr_data[obs_data$index[i]+1])
}

obs_data$viral_load <- round(pmax(pmin(obs_data$mu + rnorm(nrow(obs_data),0,sigma),45),0),6)

obs_data$viral_load <- c(0,8.2445,2.3393,18.2477,25.9738,0,12.9064,28.7029,27.1031,0,16.6276,29.0459,23.9215,25.7301,12.1522,12.6907,2.1197,18.6972,15.2154,9.7228,16.9616,39.806,22.8743,0,0,0,1.4154,20.4609,18.3101,12.7262,0,0)

obs_data$subtype <- subtype[obs_data$index+1]

log_lh_test <- 0

wp_mean_test <- c(3.5,3.3)
wp_sd_test <- c(1,0.9)
tp_sd_test <- c(1.5,1.2)
dp_mean_test <- c(36,34)
dp_sd_test <- c(3,5)
wr_mean_test <- c(5.5,5.3)
wr_sd_test <- c(2,2.1)

wp_mean_test <- c(3.5,3.3)
wp_sd_test <- c(1,1)
tp_sd_test <- c(1.5,1.5)
dp_mean_test <- c(36,34)
dp_sd_test <- c(3,3)
wr_mean_test <- c(5.5,5.3)
wr_sd_test <- c(2,2)

log_lh_test <- sum(log((sensitivity*dnorm(obs_data$viral_load-obs_data$mu,0,sd=sigma)) +
                        ((1-sensitivity)*dexp(obs_data$viral_load,1/log(10)))))

obs_data$lh_i <- log((sensitivity*dnorm(obs_data$viral_load-obs_data$mu,0,sd=sigma)) +
                       ((1-sensitivity)*dexp(obs_data$viral_load,1/log(10))))

print(log((sensitivity*dnorm(obs_data$viral_load-obs_data$mu,0,sd=sigma)) +
            ((1-sensitivity)*dexp(obs_data$viral_load,1/log(10)))))

print(log_lh_test)

for(i in 1:length(wp_data)){
  model_i <- model[i]
  index_r <- i
  subtype_i <- subtype[index_r]+1
  if(model_i == 1){
    log_lh_i <- log(dnorm(wp_data[index_r],
                          wp_mean_test[subtype_i],
                          wp_sd_test[subtype_i])) +
                log(dnorm(tp_data[index_r],
                          0,
                          tp_sd_test[subtype_i])) +
                log(dnorm(dp_data[index_r],
                          dp_mean_test[subtype_i],
                          dp_sd_test[subtype_i]))
    
    print(paste0("index ",index_r," subtype_i ",subtype_i," wp ",wp_data[index_r],
                 " wp_lh_i ",log(dnorm(wp_data[index_r],
                                                                             wp_mean_test[subtype_i],
                                                                             wp_sd_test[subtype_i])),
                " tp ", tp_data[index_r],
                 " tp_lh_i ",log(dnorm(tp_data[index_r],
                                   0,
                                   tp_sd_test[subtype_i])),
                " dp mean ",dp_mean_test[subtype_i],
                " dp sd ",dp_sd_test[subtype_i],
                " dp ", dp_data[index_r],
                " dp_lh_i ",log(dnorm(dp_data[index_r],
                                   dp_mean_test[subtype_i],
                                   dp_sd_test[subtype_i]))))
    
    log_lh_test <- log_lh_test + log_lh_i
  } else if(model_i == 2){
    log_lh_i <- log(dnorm(wp_data[index_r],
                          wp_mean_test[subtype_i],
                          wp_sd_test[subtype_i])) +
                log(dnorm(tp_data[index_r],
                          0,
                          tp_sd_test[subtype_i])) +
                log(dnorm(dp_data[index_r],
                          dp_mean_test[subtype_i],
                          dp_sd_test[subtype_i])) +
                log(dnorm(wr_data[index_r],
                          wr_mean_test[subtype_i],
                          wr_sd_test[subtype_i]))
    
    print(paste0("index ",index_r," subtype_i ",subtype_i,
                 " wp ",wp_data[index_r],
                 " wp_lh_i ",log(dnorm(wp_data[index_r],
                                     wp_mean_test[subtype_i],
                                     wp_sd_test[subtype_i])),
                 " tp ",tp_data[index_r],
                  " tp_lh_i ",log(dnorm(tp_data[index_r],
                                     0,
                                     tp_sd_test[subtype_i])),
                 " dp mean ",dp_mean_test[subtype_i],
                 " dp sd ",dp_sd_test[subtype_i],
                 " dp ",dp_data[index_r],
                  " dp_lh_i ",log(dnorm(dp_data[index_r],
                                     dp_mean_test[subtype_i],
                                     dp_sd_test[subtype_i])),
                 " wr ",wr_data[index_r],
                  " wr_lh_i ",log(dnorm(wr_data[index_r],
                                     wr_mean_test[subtype_i],
                                     wr_sd_test[subtype_i]))))
    
    log_lh_test <- log_lh_test + log_lh_i
    
  } else if(model_i == 3){
    log_lh_i <- log(dnorm(tp_data[index_r],
                          0,
                          tp_sd_test[subtype_i])) +
                log(dnorm(dp_data[index_r],
                          dp_mean_test[subtype_i],
                          dp_sd_test[subtype_i])) +
                log(dnorm(wr_data[index_r],
                          wr_mean_test[subtype_i],
                          wr_sd_test[subtype_i]))
    
    print(paste0("index ",index_r," subtype_i ",subtype_i,
                 " tp ",tp_data[index_r],
                 " tp_lh_i ",log(dnorm(tp_data[index_r],
                             0,
                             tp_sd_test[subtype_i])),
                 " dp mean ",dp_mean_test[subtype_i],
                 " dp sd ",dp_sd_test[subtype_i],
                " dp ",dp_data[index_r],
          " dp_lh_i ",log(dnorm(dp_data[index_r],
                             dp_mean_test[subtype_i],
                             dp_sd_test[subtype_i])),
          " wr ",wr_data[index_r],
          " wr_lh_i ",log(dnorm(wr_data[index_r],
                             wr_mean_test[subtype_i],
                             wr_sd_test[subtype_i]))))
  
    log_lh_test <- log_lh_test + log_lh_i
  }
  
  if(is.na(log_lh_test)){
    break
  }
}

aggregate(lh_i ~ index, data=obs_data,sum)

print(log_lh_test)
sink(file=NULL)

sink(file="lh_cpp.txt")
test_likelihood_calc()
test_likelihood_individual_calc()
sink(file=NULL)