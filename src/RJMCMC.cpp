#include <Rcpp.h>
#include "utilities.h"
#include "structs.h"
#include "make_structs.h"
#include "distributions.h"
#include "log_likelihood.h"

using namespace Rcpp;

//[[Rcpp::export]]
void rjmcmc_r(std::string output_dir,
              std::vector<int> index_arg,
              std::vector<double> viral_load_arg,
              std::vector<double> time_arg,
              std::vector<int> subtype_arg,
              std::vector<double> t_first_positive_arg,
              std::vector<double> t_last_positive_arg,
              std::vector<double> t_first_test_arg,
              std::vector<double> t_last_test_arg,
              std::vector<int> n_positive_tests_arg,
              double lod_arg,
              double sensitivity_arg,
              int n_iterations_arg,
              int n_subtypes_arg,
              int n_subjects_arg,
              int n_data_arg,
              std::vector<double> wp_current_arg,
              std::vector<double> tp_current_arg,
              std::vector<double> dp_current_arg,
              std::vector<double> wr_current_arg,
              std::vector<int> model_current_arg,
              double sigma_arg,
              std::vector<double> wp_mean_arg,
              std::vector<double> wp_sd_arg,
              std::vector<double> tp_sd_arg,
              std::vector<double> dp_mean_arg,
              std::vector<double> dp_sd_arg,
              std::vector<double> wr_mean_arg,
              std::vector<double> wr_sd_arg,
              double wp_max_arg,
              double wr_max_arg,
              double wpmean_max_arg,
              double dpmean_max_arg,
              double wrmean_max_arg,
              double wpsd_max_arg,
              double tpsd_max_arg,
              double dpsd_max_arg,
              double wrsd_max_arg,
              double sigma_max_arg,
              double wpsd_min_arg,
              double tpsd_min_arg,
              double dpsd_min_arg,
              double wrsd_min_arg,
              double sigma_min_arg,
              double p_model_1_arg,
              double p_model_2_arg,
              double p_model_3_arg,
              double wpmean_mean_arg,
              double wpmean_sd_arg,
              double dpmean_mean_arg,
              double dpmean_sd_arg,
              double wrmean_mean_arg,
              double wrmean_sd_arg,
              double wpsd_scale_arg,
              double tpsd_scale_arg,
              double dpsd_scale_arg,
              double wrsd_scale_arg,
              double sigma_scale_arg,
              std::vector<double> wp_mean_sf_arg,
              std::vector<double> wp_sd_sf_arg,
              std::vector<double> tp_sd_sf_arg,
              std::vector<double> dp_mean_sf_arg,
              std::vector<double> dp_sd_sf_arg,
              std::vector<double> wr_mean_sf_arg,
              std::vector<double> wr_sd_sf_arg,
              double sigma_sf_arg,
              double wp_sf_arg,
              double tp_sf_arg,
              double dp_sf_arg,
              double wr_sf_arg){
  // Organize data into structs
  settings_struct settings;
  viral_data_struct viral_data;
  current_data_struct current_data;
  current_parameters_struct current_parameters;
  priors_struct priors;
  scaling_factors_struct scaling_factors;
  
  make_settings_struct(settings,
                       lod_arg,
                       sensitivity_arg,
                       n_iterations_arg,
                       n_subtypes_arg,
                       n_subjects_arg,
                       n_data_arg);
  
  make_viral_data_struct(viral_data,
                         index_arg,
                         viral_load_arg, 
                         time_arg, 
                         subtype_arg,
                         t_first_positive_arg,
                         t_last_positive_arg,
                         t_first_test_arg,
                         t_last_test_arg,
                         n_positive_tests_arg);
  
  make_current_data_struct(current_data,
                           wp_current_arg,
                           tp_current_arg,
                           dp_current_arg,
                           wr_current_arg,
                           model_current_arg);
  
  make_current_parameters_struct(current_parameters,
                                 0.0,
                                 sigma_arg,
                                 wp_mean_arg,
                                 wp_sd_arg,
                                 tp_sd_arg,
                                 dp_mean_arg,
                                 dp_sd_arg,
                                 wr_mean_arg,
                                 wr_sd_arg);
  
  make_priors_struct(priors,
                     wp_max_arg,
                     wr_max_arg,
                     wpmean_max_arg,
                     dpmean_max_arg,
                     wrmean_max_arg,
                     wpsd_max_arg,
                     tpsd_max_arg,
                     dpsd_max_arg,
                     wrsd_max_arg,
                     sigma_max_arg,
                     wpsd_min_arg,
                     tpsd_min_arg,
                     dpsd_min_arg,
                     wrsd_min_arg,
                     sigma_min_arg,
                     p_model_1_arg,
                     p_model_2_arg,
                     p_model_3_arg,
                     wpmean_mean_arg,
                     wpmean_sd_arg,
                     dpmean_mean_arg,
                     dpmean_sd_arg,
                     wrmean_mean_arg,
                     wrmean_sd_arg,
                     wpsd_scale_arg,
                     tpsd_scale_arg,
                     dpsd_scale_arg,
                     wrsd_scale_arg,
                     sigma_scale_arg);
  
  make_scaling_factors_struct(scaling_factors,
                              wp_mean_sf_arg,
                              wp_sd_sf_arg,
                              tp_sd_sf_arg,
                              dp_mean_sf_arg,
                              dp_sd_sf_arg,
                              wr_mean_sf_arg,
                              wr_sd_sf_arg,
                              sigma_sf_arg,
                              wp_sf_arg,
                              tp_sf_arg,
                              dp_sf_arg,
                              wr_sf_arg);
  
  // Initialize likelihood
  current_parameters.log_likelihood = log_likelihood(viral_data,
                                                     current_data,
                                                     current_parameters,
                                                     settings);
  
  // Begin loop
  
  for(int iter = 0; iter < settings.n_iterations; iter++){
    // At each iteration:
    // Update wp_mean for each type
    // Update wp_sd for each type
    // Update dp_mean for each type
    // Update dp_sd for each type
    // Update tp_sd for each type
    // Update wr_mean for each type
    // Update wr_sd for each type
    // Update sigma
    // Update model for each individual
    // Update wp for each individual
    // Update dp for each individual
    // Update tp for each individual
    // Update wr for each individual
    // Output results
  }
  
  
}