#!/usr/bin/env Rscript
seeds <- 151:250

dir_names <- paste0("./analysis_seed_",seeds,"/")

seeds_keep <- c()

if(dir.exists(dir_names[1])){
  load(paste0(dir_names[1],"quantile_coverage.RData"))
  seeds_keep <- c(seeds_keep,seeds[1])
  loop_begin <- 2
} else{
  load(paste0(dir_names[2],"quantile_coverage.RData"))
  seeds_keep <- c(seeds_keep,seeds[2])
  loop_begin <- 3
}

coverage <- quantile_coverage

for(i in loop_begin:length(dir_names)){
  if(!dir.exists(dir_names[i])){
    next
  }
  
  seeds_keep <- c(seeds_keep,seeds[i])
  
  
  load(paste0(dir_names[i],"quantile_coverage.RData"))
  
  coverage <- rbind(coverage,quantile_coverage)
}

coverage$seed <- seeds_keep
coverage$wp_mean_cover <- coverage$wp_mean_0_cover + coverage$wp_mean_1_cover + coverage$wp_mean_2_cover
coverage$dp_mean_cover <- coverage$dp_mean_0_cover + coverage$dp_mean_1_cover + coverage$dp_mean_2_cover
coverage$wr_mean_cover <- coverage$wr_mean_0_cover + coverage$wr_mean_1_cover + coverage$wr_mean_2_cover
coverage$n_cover <- coverage$wp_mean_cover + coverage$dp_mean_cover + coverage$wr_mean_cover

coverage$keep <- coverage$wp_mean_cover > 1 & coverage$dp_mean_cover > 1 & coverage$wr_mean_cover > 1

coverage$final_keep <- (coverage$n_cover > 7 & coverage$keep) | (coverage$n_cover == 7 & coverage$dp_sd_0_cover) | (coverage$n_cover == 7 & coverage$dp_mean_0_cover)

save(coverage,file=paste0("coverage_seeds_",min(seeds_keep),"_",max(seeds_keep),".RData"))
