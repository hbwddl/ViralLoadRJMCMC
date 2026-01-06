library(dplyr)

load("~/Documents/Research/Within-Host/RJMCMC_Results/multi_chain_flu_summary_seeds_1111_1115/dp_out.RData")
load("~/Documents/Research/Within-Host/RJMCMC_Results/multi_chain_flu_summary_seeds_1111_1115/model_out.RData")
load("~/Documents/Research/Within-Host/RJMCMC_Results/multi_chain_flu_summary_seeds_1111_1115/scalars_out.RData")
load("~/Documents/Research/Within-Host/RJMCMC_Results/multi_chain_flu_summary_seeds_1111_1115/tp_out.RData")
load("~/Documents/Research/Within-Host/RJMCMC_Results/multi_chain_flu_summary_seeds_1111_1115/wp_out.RData")
load("~/Documents/Research/Within-Host/RJMCMC_Results/multi_chain_flu_summary_seeds_1111_1115/wr_out.RData")

setwd("~/Documents/Research/Within-Host/RJMCMC_Results/analysis_flu_seed_1111")

load("individual_data_in.RData")
load("priors.RData")
load("settings_in.RData")
load("viral_data_in.RData")

setwd("~/Documents/Research/Within-Host/RJMCMC_Results/multi_chain_flu_summary_seeds_1111_1115")

lod = 45

n_graph_rows = 5
n_graph_cols = 4

graph_img_starts = seq(from=1,to=nrow(individual_data),by=n_graph_cols*n_graph_rows)
graph_img_ends = sort(unique(c(seq(from=n_graph_rows*n_graph_cols,to=nrow(individual_data),by=n_graph_cols*n_graph_rows),nrow(individual_data))))

for(j in 1:length(graph_img_starts)){
  png(filename=paste0("data_plots_",j,".png"),width = 2400, height = 3000,res=200)
  
  par(mfrow=c(n_graph_rows,n_graph_cols))
  
  for(i in (graph_img_starts[j]):(graph_img_ends[j])){
    plot_dat <- viral_data %>%
      filter(index==i-1) %>%
      arrange(time)
    
    model_infer <- which.max(table(c(model_out[,i],1,2,3)))
    
    which_infer_model <- model_out[,i] == model_infer
    
    lower_bound_adj = individual_data$t_first_test[i]
    last_day_fair = individual_data$t_last_test[i] - lower_bound_adj
    
    wp_quantile <- quantile(wp_out[which_infer_model,i],probs = c(0.025,0.5,0.975))
    tp_quantile <- quantile(tp_out[which_infer_model,i],probs = c(0.025,0.5,0.975))
    dp_quantile <- quantile(dp_out[which_infer_model,i],probs = c(0.025,0.5,0.975))
    wr_quantile <- quantile(wr_out[which_infer_model,i],probs = c(0.025,0.5,0.975))
    
    tp_med <- median(tp_out[which_infer_model,i])
    
    plot_col <- c("blue","red","purple")[individual_data$subtype[i]+1]
    
    plot(plot_dat$time - lower_bound_adj,lod - plot_dat$viral_load,pch=19,main=paste0("Observed Data, ID ",individual_data$pig_id[i]),col=plot_col,
         xlim=c(min(-wp_quantile[3] - lower_bound_adj,plot_dat$time - lower_bound_adj),max(wr_quantile[3],plot_dat$time)),
         # ylim=c(0,max(dp_quantile[3],plot_dat$time)),
         ylim=c(lod,lod-max(dp_quantile[3],plot_dat$time)),
         xlab="Time",ylab="Ct")
    
    if(model_infer != 3){
      lines(c(tp_med-wp_quantile[2] - lower_bound_adj,tp_med - lower_bound_adj),c(lod,lod-dp_quantile[2]))
      lines(c(tp_med-wp_quantile[1] - lower_bound_adj,tp_med - lower_bound_adj),c(lod,lod-dp_quantile[1]),lty="dashed")
      lines(c(tp_med-wp_quantile[3] - lower_bound_adj,tp_med - lower_bound_adj),c(lod,lod-dp_quantile[3]),lty="dashed")
    }
    
    if(model_infer != 1){
      lines(c(tp_med - lower_bound_adj,tp_med+wr_quantile[2] - lower_bound_adj),c(lod-dp_quantile[2],lod))
      lines(c(tp_med - lower_bound_adj,tp_med+wr_quantile[1] - lower_bound_adj),c(lod-dp_quantile[1],lod),lty="dashed")
      lines(c(tp_med - lower_bound_adj,tp_med+wr_quantile[3] - lower_bound_adj),c(lod-dp_quantile[3],lod),lty="dashed")
    }
    
    abline(v = last_day_fair, lty = "dashed", col = "grey80")
  }
  
  par(mfrow=c(1,1))
  
  dev.off()
}

