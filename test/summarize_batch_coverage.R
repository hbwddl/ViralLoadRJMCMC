#!/usr/bin/env Rscript
seeds <- 1:50

dir_names <- paste0("./analysis_seed_",seeds,"/")

load(paste0(dir_names[1],"quantile_coverage.RData"))

coverage <- quantile_coverage

for(i in 2:length(dir_names)){
  load(paste0(dir_names[i],"quantile_coverage.RData"))
  
  coverage <- rbind(coverage,quantile_coverage)
}

save(coverage,file=paste0("coverage_seeds_",min(seeds),"_",max(seeds)))