load("~/Documents/Research/Within-Host/RJMCMC_Results/analysis_flu_seed_1119/individual_data_in.RData")
load("~/Documents/Research/Within-Host/RJMCMC_Results/analysis_flu_seed_1119/viral_data_in.RData")

n_gt0 <- aggregate(viral_data$viral_load ~ viral_data$index,data=viral_data,function(x){return(sum(x>0))})

individual_data$n_gt0 <- n_gt0[match(n_gt0[,1],individual_data$index),2]

aggregate(individual_data$n_gt0 ~ individual_data$subtype,data=individual_data,mean)

n_test <- aggregate(viral_data$viral_load ~ viral_data$index,data=viral_data,length)

individual_data$n_test <- n_test[match(n_test[,1],individual_data$index),2]

aggregate(individual_data$n_test ~ individual_data$subtype,data=individual_data,mean)

first_gt0 <- aggregate(viral_data$viral_load ~ viral_data$index,data=viral_data,function(x){return(min(which(x>0)))})

individual_data$first_gt0 <- first_gt0[match(first_gt0[,1],individual_data$index),2]

aggregate(individual_data$first_gt0 ~ individual_data$subtype,data=individual_data,mean)
