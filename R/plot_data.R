plot_sim_data <- function(sim_data_arg,
                          outfile_name="./sim_plots.pdf"){
  plot_xlim <- c(min(sim_data_arg$individual_data$tp-sim_data_arg$individual_data$wp),max(sim_data_arg$individual_data$tp+sim_data_arg$individual_data$wr))
  plot_ylim <- c(0,sim_data_arg$settings$lod)
  
  pdf(file=outfile_name,width=10,height=10)
  par(mfrow=c(3,3))
  for(i in 1:nrow(sim_data_arg$individual_data)){
    plot_data_viral <- sim_data_arg$viral_load_data %>%
                        filter(index_r == i)
    
    plot(plot_data_viral$time,
         plot_data_viral$viral_load_obs,
         xlim=plot_xlim,
         ylim=plot_ylim,main=paste0("Index ",i),
         pch=19)
    abline(v=sim_data_arg$individual_data$t_first_test[i],lty=2)
    abline(v=sim_data_arg$individual_data$t_last_test[i],lty=2)
    segments(sim_data_arg$individual_data$tp[i]-sim_data_arg$individual_data$wp[i],
             0,
             sim_data_arg$individual_data$tp[i],
             sim_data_arg$individual_data$dp[i])
    segments(sim_data_arg$individual_data$tp[i],
             sim_data_arg$individual_data$dp[i],
             sim_data_arg$individual_data$tp[i]+sim_data_arg$individual_data$wr[i],
             0)
    
  }
  par(mfrow=c(1,1))
  dev.off()
}