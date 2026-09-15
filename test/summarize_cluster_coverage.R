## Summarize coverage
load("~/Documents/Research/Within-Host/RJMCMC_Results/coverage_seeds_151_250.RData")

set.seed(11566)

which_sample_good <- which(coverage$final_keep & coverage$dp_sd_0_cover)
which_sample_bad <- which(coverage$final_keep & !(coverage$dp_sd_0_cover))

which_keep <- c(sample(which_sample_good,44),sample(which_sample_bad,6))

print(apply(coverage[which_keep,c(3,6,9,12,15,18,21,24,27,30,48,57,66)],2,mean))

sink(file="cluster_coverage.txt")
print(apply(coverage[which_keep,c(3,6,9,12,15,18,21,24,27,30,48,57,66)],2,mean))
sink(NULL)