ntest <- 1000

rnorm_results <- rep(NA,ntest)

for(i in 1:ntest){
  rnorm_results[i] <- test_rnorm_boost(i)
}

hist(rnorm_results)

runif_results <- rep(NA,ntest)

for(i in 1:ntest){
  runif_results[i] <- test_runif(i)
}

hist(runif_results)
