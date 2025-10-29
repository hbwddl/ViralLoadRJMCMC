## R RJMCMC function
log_likelihood_r <- function(viral_df,
                              individual_df,
                              settings_vec,
                              wp_mean_vec,
                              wp_sd_vec,
                              dp_mean_vec,
                              dp_sd_vec,
                              tp_sd_vec,
                              wr_mean_vec,
                              wr_sd_vec,
                              wp_vec,
                              dp_vec,
                              tp_vec,
                              wr_vec,
                              model_vec,
                              sigma){
  
  index_vec <- viral_df[,"index"]
  viral_load_vec <- viral_df[,"viral_load"]
  time_vec <- viral_df[,"time"]
  subtype_vec <- individual_df[,"subtype"]
  t_first_positive_vec <- individual_df[,"t_first_positive"]
  t_last_positive_vec <- individual_df[,"t_last_positive"]
  t_first_test_vec <- individual_df[,"t_first_test"]
  t_last_test_vec <- individual_df[,"t_last_test"]
  n_positive_tests <- individual_df[,"n_positive"]
  max_viral_load <- individual_df[,"max_viral_load"]

  lod_in <- settings_vec[,"lod"]
  sensitivity_in <- settings_vec[,"sensitivity"]
  n_subtypes_in <- settings_vec[,"n_subtypes"]
  n_subjects_in <- settings_vec[,"n_subjects"]
  n_data_in <- settings_vec[,"n_data"]
  
  
  
  return(log_likelihood_r_exp(index_vec,
           viral_load_vec,
           time_vec,
           subtype_vec,
           t_first_positive_vec,
           t_last_positive_vec,
           t_first_test_vec,
           t_last_test_vec,
           n_positive_tests,
           max_viral_load,
           lod_in,
           sensitivity_in,
           n_subtypes_in,
           n_subjects_in,
           n_data_in,
           wp_vec,
           tp_vec,
           dp_vec,
           wr_vec,
           model_vec,
           sigma,
           wp_mean_vec,
           wp_sd_vec,
           tp_sd_vec,
           dp_mean_vec,
           dp_sd_vec,
           wr_mean_vec,
           wr_sd_vec
           ));
}