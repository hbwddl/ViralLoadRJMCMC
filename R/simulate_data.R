## Simulating data
simulate_viral_load_data <- function(data_settings_arg,
                                     model_parameters_arg,
                                     sim_seed=as.integer(Sys.time())){
  require(dplyr)
  set.seed(sim_seed)
  
  lod <- data_settings_arg[["lod"]]
  n_pop <- data_settings_arg[["n"]]
  p_group <- data_settings_arg[["p_group"]]
  n_group <- length(p_group)
  t_obs <- data_settings_arg[["t_obs"]]
  sensitivity <- data_settings_arg[["sensitivity"]]
  
  p_model <- param_settings_in[["p_model"]]
  wp_mean <- param_settings_in[["wp_mean"]]
  wp_sd <- param_settings_in[["wp_sd"]]
  tp_sd <- param_settings_in[["tp_sd"]]
  dp_mean <- param_settings_in[["dp_mean"]]
  dp_sd <- param_settings_in[["dp_sd"]]
  wr_mean <- param_settings_in[["wr_mean"]]
  wr_sd <- param_settings_in[["wr_sd"]]
  sigma_val <- param_settings_in[["sigma"]]
  wp_min <- param_settings_in[["wp_min"]]
  wp_max <- param_settings_in[["wp_max"]]
  dp_min <- param_settings_in[["dp_min"]]
  tp_min <- param_settings_in[["tp_min"]]
  tp_max <- param_settings_in[["tp_max"]]
  wr_min <- param_settings_in[["wr_min"]]
  wr_max <- param_settings_in[["wr_max"]]
  wp_mean <- param_settings_in[["wp_mean"]]
  
  indiv_data <- data.frame(index=0:(n_pop-1),
                           index_r=1:n_pop,
                           group_r=sample(1:length(p_group),n_pop,replace=T,prob=p_group),
                           model=sample(1:3,n_pop,replace=T,prob=p_model))
  
  indiv_data$group <- indiv_data$group_r-1
  
  indiv_data$wp <- truncnorm::rtruncnorm(n_pop,a=wp_min,b=wp_max,mean=wp_mean[indiv_data$group_r],sd=wp_sd[indiv_data$group_r])
  indiv_data$tp_actual <- rnorm(n_pop,0,tp_sd[indiv_data$group_r])
  indiv_data$dp <- truncnorm::rtruncnorm(n_pop,a=dp_min,b=lod,mean=dp_mean[indiv_data$group_r],sd=dp_sd[indiv_data$group_r])
  indiv_data$wr <- truncnorm::rtruncnorm(n_pop,a=wr_min,b=wr_max,mean=wr_mean[indiv_data$group_r],sd=wr_sd[indiv_data$group_r])
  
  # indiv_data$tp_actual <- ifelse(indiv_data$model==1,
  #                         runif(1,min=max(t_obs),max=max(t_obs)+(wp_max/2)),
  #                         ifelse(indiv_data$model==3,
  #                                runif(1,min=min(t_obs)-(wr_max/2),max=min(t_obs)),
  #                                runif(1,min=min(t_obs),max=max(t_obs))))
  
  # for(i in 1:nrow(indiv_data)){
  #   indiv_data$tp_actual[i] <- ifelse(indiv_data$model[i]==1,
  #                              runif(1,min=max(t_obs),max=max(t_obs)+(wp_max*0.75)),
  #                              ifelse(indiv_data$model[i]==3,
  #                                     runif(1,min=min(t_obs)-(wr_max*0.75),max=min(t_obs)),
  #                                     runif(1,min=min(t_obs),max=max(t_obs))))
  #   
  #   if(indiv_data$model[i] == 1 & indiv_data$tp_actual[i] - indiv_data$wp[i] > (max(t_obs))){
  #     indiv_data$tp_actual[i] <- runif(1,min=max(t_obs),max=max(t_obs)+(indiv_data$wp[i])/20)
  #   }
  #   if(indiv_data$model[i] == 3 & indiv_data$tp_actual[i] + indiv_data$wr[i] < (min(t_obs))){
  #     indiv_data$tp_actual[i] <- runif(1,min=min(t_obs) - (indiv_data$wr[i]/20),max=min(t_obs))
  #   }
  # }
  
  viral_data_g <- expand.grid(indiv_data$index_r,
                              t_obs)
  
  viral_data <- data.frame(index_r = viral_data_g[,1],
                           time_actual = viral_data_g[,2])
  
  viral_data$index <- viral_data$index_r - 1
  viral_data$viral_load_obs <- 0
  
  for(i in 1:nrow(viral_data)){
    viral_data$viral_load_mu[i] <- mu(viral_data$time_actual[i], 
                                   indiv_data$wp[viral_data$index_r[i]], 
                                   indiv_data$tp_actual[viral_data$index_r[i]], 
                                   indiv_data$dp[viral_data$index_r[i]], 
                                   indiv_data$wr[viral_data$index_r[i]])
  }
  
  viral_data$viral_load_obs <- NA
  
  for(i in 1:nrow(viral_data)){
    if(viral_data$viral_load_mu[i] == 0){
      runif_false_pos <- runif(1,0,1)
      if(runif_false_pos > sensitivity){
        viral_data$viral_load_obs[i] <- rexp(1,1/log(10))
      } else{
        viral_data$viral_load_obs[i] <- 0
      }
    } else{
      viral_data$viral_load_obs[i] <- min(max(viral_data$viral_load_mu[i] + rnorm(1,0,sigma_val),0),lod-2)
    }
  }
    
  indiv_data$day_peak_obs <- 0
  indiv_data$first_test_actual <- 0
  indiv_data$last_test_actual <- 0
  indiv_data$first_gt0_actual <- 0
  indiv_data$last_gt0_actual <- 0
  indiv_data$tp <- 0
  indiv_data$n_positive <- 0
  indiv_data$max_viral_load <- 0
  
  viral_model <- indiv_data$model[viral_data$index_r]
  viral_tp_actual <- indiv_data$tp_actual[viral_data$index_r]
  
  viral_keep <- T
  
  for(i in 1:nrow(viral_data)){
    viral_keep[i] <- ifelse(viral_model[i] == 2,
                                        T,
                                        ifelse(viral_model[i] == 1,
                                               viral_data$time_actual[i] <= viral_tp_actual[i],
                                               ifelse(viral_model[i] == 3,
                                               viral_data$time_actual[i] >= viral_tp_actual[i],
                                                UNTITLED())))
  }
  
  viral_data <- viral_data %>%
                  filter(viral_keep)
  
  indiv_data <- indiv_data %>%
                  filter(index_r %in% viral_data$index_r)
  
  indiv_index_r_adj <- 1:nrow(indiv_data)
  viral_index_r_adj <- indiv_index_r_adj[match(viral_data$index_r,indiv_data$index_r)]

  # indiv_data$indiv_index_r_adj <- indiv_index_r_adj
  # viral_data$viral_index_r_adj <- viral_index_r_adj
  
  indiv_index_adj <- indiv_index_r_adj - 1
  viral_index_adj <- viral_index_r_adj - 1

  indiv_data$index <- indiv_index_adj
  viral_data$index <- viral_index_adj

  indiv_data$index_r <- indiv_index_r_adj
  viral_data$index_r <- viral_index_r_adj

  ## Adjust time to 0 at peak value
  for(i in 1:nrow(indiv_data)){
    viral_data_i <- viral_data %>%
                    filter(index_r == indiv_data$index_r[i]) %>%
                    arrange(time_actual)
    
      peak_day <- viral_data_i$time_actual[which.max(viral_data_i$viral_load_obs)]
      first_test <- min(viral_data_i$time_actual)
      last_test <- max(viral_data_i$time_actual)
      first_gt0 <- viral_data_i$time_actual[min(which(viral_data_i$viral_load_obs > 0))]
      last_gt0 <- viral_data_i$time_actual[max(which(viral_data_i$viral_load_obs > 0))]
      max_viral_load <- max(viral_data_i$viral_load_obs)
      
      indiv_data$day_peak_obs[i] <- peak_day
      indiv_data$first_test_actual[i] <- first_test
      indiv_data$last_test_actual[i] <- last_test
      indiv_data$first_gt0_actual[i] <- first_gt0
      indiv_data$last_gt0_actual[i] <- last_gt0
      indiv_data$n_positive[i] <- sum(viral_data_i$viral_load_obs > 0)
      indiv_data$max_viral_load[i] <- max_viral_load
    
  }
  
  indiv_data$tp <- indiv_data$tp_actual - indiv_data$day_peak_obs
  indiv_data$t_first_positive <- indiv_data$first_gt0_actual - indiv_data$day_peak_obs
  indiv_data$t_last_positive <- indiv_data$last_gt0_actual - indiv_data$day_peak_obs
  indiv_data$t_first_test <- indiv_data$first_test_actual - indiv_data$day_peak_obs
  indiv_data$t_last_test <- indiv_data$last_test_actual - indiv_data$day_peak_obs
  
  viral_data$time <- NA
  
  for(i in 1:nrow(viral_data)){
    viral_data$time[i] <- viral_data$time_actual[i] - indiv_data$day_peak_obs[viral_data$index_r[i]]
  }
  
  sim_data <- list(settings = data_settings_arg,
                   parameters = model_parameters_arg,
                   individual_data = indiv_data,
                   viral_load_data = viral_data,
                   seed = mcmc_seed)
  
  return(sim_data)
}

## Simulating data
simple_simulate_viral_load_data <- function(data_settings_arg,
                                     model_parameters_arg,
                                     sim_seed=as.integer(Sys.time())){
  require(dplyr)
  set.seed(sim_seed)
  
  lod <- data_settings_arg[["lod"]]
  n_pop <- data_settings_arg[["n"]]
  p_group <- data_settings_arg[["p_group"]]
  n_group <- length(p_group)
  t_obs <- data_settings_arg[["t_obs"]]
  sensitivity <- data_settings_arg[["sensitivity"]]
  
  p_model <- param_settings_in[["p_model"]]
  wp_mean <- param_settings_in[["wp_mean"]]
  wp_sd <- param_settings_in[["wp_sd"]]
  tp_sd <- param_settings_in[["tp_sd"]]
  dp_mean <- param_settings_in[["dp_mean"]]
  dp_sd <- param_settings_in[["dp_sd"]]
  wr_mean <- param_settings_in[["wr_mean"]]
  wr_sd <- param_settings_in[["wr_sd"]]
  sigma_val <- param_settings_in[["sigma"]]
  wp_min <- param_settings_in[["wp_min"]]
  wp_max <- param_settings_in[["wp_max"]]
  dp_min <- param_settings_in[["dp_min"]]
  tp_min <- param_settings_in[["tp_min"]]
  tp_max <- param_settings_in[["tp_max"]]
  wr_min <- param_settings_in[["wr_min"]]
  wr_max <- param_settings_in[["wr_max"]]
  wp_mean <- param_settings_in[["wp_mean"]]
  
  indiv_data <- data.frame(index=0:(n_pop-1),
                           index_r=1:n_pop,
                           group_r=sample(1:length(p_group),n_pop,replace=T,prob=p_group),
                           model=sample(1:3,n_pop,replace=T,prob=p_model))
  
  indiv_data$group <- indiv_data$group_r-1
  
  indiv_data$wp <- truncnorm::rtruncnorm(n_pop,a=wp_min,b=wp_max,mean=wp_mean[indiv_data$group_r],sd=wp_sd[indiv_data$group_r])
  indiv_data$tp <- rnorm(n_pop,0,tp_sd[indiv_data$group_r])
  indiv_data$dp <- truncnorm::rtruncnorm(n_pop,a=dp_min,b=lod,mean=dp_mean[indiv_data$group_r],sd=dp_sd[indiv_data$group_r])
  indiv_data$wr <- truncnorm::rtruncnorm(n_pop,a=wr_min,b=wr_max,mean=wr_mean[indiv_data$group_r],sd=wr_sd[indiv_data$group_r])
  
  viral_data_g <- expand.grid(indiv_data$index_r,
                              t_obs)
  
  viral_data <- data.frame(index_r = viral_data_g[,1],
                           time_actual = viral_data_g[,2])
  
  viral_data$index <- viral_data$index_r - 1
  viral_data$viral_load_obs <- 0
  
  for(i in 1:nrow(viral_data)){
    viral_data$viral_load_mu[i] <- mu(viral_data$time_actual[i], 
                                      indiv_data$wp[viral_data$index_r[i]], 
                                      indiv_data$tp_actual[viral_data$index_r[i]], 
                                      indiv_data$dp[viral_data$index_r[i]], 
                                      indiv_data$wr[viral_data$index_r[i]])
  }
  
  viral_data$viral_load_obs <- NA
  
  for(i in 1:nrow(viral_data)){
    if(viral_data$viral_load_mu[i] == 0){
      runif_false_pos <- runif(1,0,1)
      if(runif_false_pos > sensitivity){
        viral_data$viral_load_obs[i] <- rexp(1,1/log(10))
      } else{
        viral_data$viral_load_obs[i] <- 0
      }
    } else{
      viral_data$viral_load_obs[i] <- min(max(viral_data$viral_load_mu[i] + rnorm(1,0,sigma_val),0),lod-2)
    }
  }
  
  indiv_data$day_peak_obs <- 0
  indiv_data$first_test_actual <- 0
  indiv_data$last_test_actual <- 0
  indiv_data$first_gt0_actual <- 0
  indiv_data$last_gt0_actual <- 0
  indiv_data$tp <- 0
  indiv_data$n_positive <- 0
  indiv_data$max_viral_load <- 0
  
  ## Adjust time to 0 at peak value
  for(i in 1:nrow(indiv_data)){
    viral_data_i <- viral_data %>%
      filter(index_r == indiv_data$index_r[i]) %>%
      arrange(time_actual)
    
    peak_day <- viral_data_i$time_actual[which.max(viral_data_i$viral_load_obs)]
    first_test <- min(viral_data_i$time_actual)
    last_test <- max(viral_data_i$time_actual)
    first_gt0 <- viral_data_i$time_actual[min(which(viral_data_i$viral_load_obs > 0))]
    last_gt0 <- viral_data_i$time_actual[max(which(viral_data_i$viral_load_obs > 0))]
    max_viral_load <- max(viral_data_i$viral_load_obs)
    
    indiv_data$day_peak_obs[i] <- peak_day
    indiv_data$first_test_actual[i] <- first_test
    indiv_data$last_test_actual[i] <- last_test
    indiv_data$first_gt0_actual[i] <- first_gt0
    indiv_data$last_gt0_actual[i] <- last_gt0
    indiv_data$n_positive[i] <- sum(viral_data_i$viral_load_obs > 0)
    indiv_data$max_viral_load[i] <- max_viral_load
  }
  
  indiv_data$tp <- indiv_data$tp_actual - indiv_data$day_peak_obs
  indiv_data$t_first_positive <- indiv_data$first_gt0_actual - indiv_data$day_peak_obs
  indiv_data$t_last_positive <- indiv_data$last_gt0_actual - indiv_data$day_peak_obs
  indiv_data$t_first_test <- indiv_data$first_test_actual - indiv_data$day_peak_obs
  indiv_data$t_last_test <- indiv_data$last_test_actual - indiv_data$day_peak_obs
  viral_data$time <- viral_data$time_actual - indiv_data$day_peak_obs[viral_data$index_r]
  
  indiv_data$index <- 0:(nrow(indiv_data)-1)
  viral_data$index <- indiv_data$index[match(viral_data$index_r,indiv_data$index_r)]
  
  indiv_data$index_r <- 1:nrow(indiv_data)
  viral_data$index_r <- indiv_data$index_r[match(viral_data$index,indiv_data$index)]
  
  sim_data <- list(settings = data_settings_arg,
                   parameters = model_parameters_arg,
                   individual_data = indiv_data,
                   viral_load_data = viral_data,
                   seed = mcmc_seed)
  
  return(sim_data)
}