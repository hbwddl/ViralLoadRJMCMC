### Log likelihood for toy data
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

obs_data$subtype <- subtype[obs_data$index+1]

log_lh_test <- 0

wp_mean_test <- c(3.5,3.3)
wp_sd_test <- c(1,0.9)
tp_sd_test <- c(1.5,1.2)
dp_mean_test <- c(36,34)
dp_sd_test <- c(3,5)
wr_mean_test <- c(5.5,5.3)
wr_sd_test <- c(2,2.1)

log_lh_test <- sum(log((sensitivity*dnorm(obs_data$viral_load-obs_data$mu,0,sd=sigma)) +
                        ((1-sensitivity)*dexp(obs_data$viral_load,1/log(10)))))

print(log((sensitivity*dnorm(obs_data$viral_load-obs_data$mu,0,sd=sigma)) +
            ((1-sensitivity)*dexp(obs_data$viral_load,1/log(10)))))

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
    
    print(paste0("index ",index_r," subtype_i ",subtype_i," wp_i ",log(dnorm(wp_data[index_r],
                                                                             wp_mean_test[subtype_i],
                                                                             wp_sd_test[subtype_i])),
                " tp_i ",log(dnorm(tp_data[index_r],
                                   0,
                                   tp_sd_test[subtype_i])),
                " dp_i ",log(dnorm(dp_data[index_r],
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
    
    print(paste0("index ",index_r," subtype_i ",subtype_i," wp_i ",log(dnorm(wp_data[index_r],
                                                                             wp_mean_test[subtype_i],
                                                                             wp_sd_test[subtype_i])),
                  " tp_i ",log(dnorm(tp_data[index_r],
                                     0,
                                     tp_sd_test[subtype_i])),
                  " dp_i ",log(dnorm(dp_data[index_r],
                                     dp_mean_test[subtype_i],
                                     dp_sd_test[subtype_i])),
                  " wr_i ",log(dnorm(wr_data[index_r],
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
    
    print(paste0("index ",index_r," subtype_i ",subtype_i," tp_i ",log(dnorm(tp_data[index_r],
                             0,
                             tp_sd_test[subtype_i])),
          " dp_i ",log(dnorm(dp_data[index_r],
                             dp_mean_test[subtype_i],
                             dp_sd_test[subtype_i])),
          " wr_i ",log(dnorm(wr_data[index_r],
                             wr_mean_test[subtype_i],
                             wr_sd_test[subtype_i]))))
  
    log_lh_test <- log_lh_test + log_lh_i
  }
  
  if(is.na(log_lh_test)){
    break
  }
}

print(log_lh_test)
sink(file=NULL)

sink(file="lh_cpp.txt")
test_likelihood_calc()
sink(file=NULL)