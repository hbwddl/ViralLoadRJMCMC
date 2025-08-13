## R RJMCMC function
viral_load_rjmcmc <- function(output_dir,
                              viral_df,
                              individual_df,
                              settings_vec,
                              wp_mean_init_vec,
                              wp_sd_init_vec,
                              dp_mean_init_vec,
                              dp_sd_init_vec,
                              tp_sd_init_vec,
                              wr_mean_init_vec,
                              wr_sd_init_vec,
                              wp_init_vec,
                              dp_init_vec,
                              tp_init_vec,
                              wr_init_vec,
                              model_init_vec,
                              sigma_init,
                              priors_vec,
                              wp_mean_sf_vec,
                              wp_sd_sf_vec,
                              dp_mean_sf_vec,
                              dp_sd_sf_vec,
                              tp_sd_sf_vec,
                              wr_mean_sf_vec,
                              wr_sd_sf_vec,
                              sigma_sf,
                              wp_sf,
                              tp_sf,
                              dp_sf,
                              wr_sf,
                              seed_arg){
  
  index_vec <- viral_df[,"index"]
  viral_load_vec <- viral_df[,"viral_load"]
  time_vec <- viral_df[,"time"]
  subtype_vec <- individual_df[,"subtype"]
  t_first_positive_vec <- individual_df[,"t_first_positive"]
  t_last_positive_vec <- individual_df[,"t_last_positive"]
  t_first_test_vec <- individual_df[,"t_first_test"]
  t_last_test_vec <- individual_df[,"t_last_test"]
  n_positive_tests <- individual_df[,"n_positive"]
  
  lod_in <- settings_vec[,"lod"]
  sensitivity_in <- settings_vec[,"sensitivity"]
  n_iterations_in <- settings_vec[,"n_iterations"]
  n_subtypes_in <- settings_vec[,"n_subtypes"]
  n_subjects_in <- settings_vec[,"n_subjects"]
  n_data_in <- settings_vec[,"n_data"]
  
  wp_max_in <- priors_vec[,"wp_max"]
  wr_max_in <- priors_vec[,"wr_max"]
  wpmean_max_in <- priors_vec[,"wpmean_max"]
  dpmean_max_in <- priors_vec[,"dpmean_max"]
  wrmean_max_in <- priors_vec[,"wrmean_max"]
  wpsd_max_in <- priors_vec[,"wpsd_max"]
  tpsd_max_in <- priors_vec[,"tpsd_max"]
  dpsd_max_in <- priors_vec[,"dpsd_max"]
  wrsd_max_in <- priors_vec[,"wrsd_max"]
  sigma_max_in <- priors_vec[,"sigma_max"]
  wpsd_min_in <- priors_vec[,"wpsd_min"]
  tpsd_min_in <- priors_vec[,"tpsd_min"]
  dpsd_min_in <- priors_vec[,"dpsd_min"]
  wrsd_min_in <- priors_vec[,"wrsd_min"]
  sigma_min_in <- priors_vec[,"sigma_min"]
  p_model_1_in <- 0.334
  p_model_2_in <- 0.333
  p_model_3_in <- 0.333
  wpmean_mean_in <- priors_vec[,"wpmean_mean"]
  wpmean_sd_in <- priors_vec[,"wpmean_sd"]
  dpmean_mean_in <- priors_vec[,"dpmean_mean"]
  dpmean_sd_in <- priors_vec[,"dpmean_sd"]
  wrmean_mean_in <- priors_vec[,"wrmean_mean"]
  wrmean_sd_in <- priors_vec[,"wrmean_sd"]
  wpsd_scale_in <- priors_vec[,"wpsd_scale"]
  tpsd_scale_in <- priors_vec[,"tpsd_scale"]
  dpsd_scale_in <- priors_vec[,"dpsd_scale"]
  wrsd_scale_in <- priors_vec[,"wrsd_scale"]
  sigma_scale_in <- priors_vec[,"sigma_scale"]
  
  rjmcmc_r(output_dir,
           index_vec,
           viral_load_vec,
           time_vec,
           subtype_vec,
           t_first_positive_vec,
           t_last_positive_vec,
           t_first_test_vec,
           t_last_test_vec,
           n_positive_tests,
           lod_in,
           sensitivity_in,
           n_iterations_in,
           n_subtypes_in,
           n_subjects_in,
           n_data_in,
           wp_init_vec,
           tp_init_vec,
           dp_init_vec,
           wr_init_vec,
           model_init_vec,
           sigma_init,
           wp_mean_init_vec,
           wp_sd_init_vec,
           tp_sd_init_vec,
           dp_mean_init_vec,
           dp_sd_init_vec,
           wr_mean_init_vec,
           wr_sd_sf_vec,
           wp_max_in,
           wr_max_in,
           wpmean_max_in,
           dpmean_max_in,
           wrmean_max_in,
           wpsd_max_in,
           tpsd_max_in,
           dpsd_max_in,
           wrsd_max_in,
           sigma_max_in,
           wpsd_min_in,
           tpsd_min_in,
           dpsd_min_in,
           wrsd_min_in,
           sigma_min_in,
           p_model_1_in,
           p_model_2_in,
           p_model_3_in,
           wpmean_mean_in,
           wpmean_sd_in,
           dpmean_mean_in,
           dpmean_sd_in,
           wrmean_mean_in,
           wrmean_sd_in,
           wpsd_scale_in,
           tpsd_scale_in,
           dpsd_scale_in,
           wrsd_scale_in,
           sigma_scale_in,
           wp_mean_sf_vec,
           wp_sd_sf_vec,
           tp_sd_sf_vec,
           dp_mean_sf_vec,
           dp_sd_sf_vec,
           wr_mean_sf_vec,
           wr_sd_sf_vec,
           sigma_sf,
           wp_sf,
           tp_sf,
           dp_sf,
           wr_sf,
           seed_arg);
}