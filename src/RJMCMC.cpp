#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include "utilities.h"
#include "structs.h"
#include "make_structs.h"
#include "distributions.h"
#include "log_likelihood.h"
#include "update_scalars.h"
#include "debug.h"

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
              std::vector<double> priors_vec,
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
              double wr_sf_arg,
              double seed_arg){
  // Organize data into structs
  settings_struct settings;
  viral_data_struct viral_data;
  current_data_struct current_data;
  current_parameters_struct current_parameters;
  priors_struct priors;
  scaling_factors_struct scaling_factors;
  
  print_pos(__FILE__,__LINE__);
  
  rng_type rng_value(seed_arg);
  
  print_pos(__FILE__,__LINE__);
  
  make_settings_struct(settings,
                       lod_arg,
                       sensitivity_arg,
                       n_iterations_arg,
                       n_subtypes_arg,
                       n_subjects_arg,
                       n_data_arg);
  
  print_pos(__FILE__,__LINE__);
  
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
  
  print_pos(__FILE__,__LINE__);
  
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
  
  print_pos(__FILE__,__LINE__);
  
  double wp_max_arg = priors_vec.at(0);
  double wr_max_arg = priors_vec.at(1);
  double wpmean_max_arg = priors_vec.at(2);
  double dpmean_max_arg = priors_vec.at(3);
  double wrmean_max_arg = priors_vec.at(4);
  double wpsd_max_arg = priors_vec.at(5);
  double tpsd_max_arg = priors_vec.at(6);
  double dpsd_max_arg = priors_vec.at(7);
  double wrsd_max_arg = priors_vec.at(8);
  double sigma_max_arg = priors_vec.at(9);
  double wpsd_min_arg = priors_vec.at(10);
  double tpsd_min_arg = priors_vec.at(11);
  double dpsd_min_arg = priors_vec.at(12);
  double wrsd_min_arg = priors_vec.at(13);
  double sigma_min_arg = priors_vec.at(14);
  double p_model_1_arg = priors_vec.at(15);
  double p_model_2_arg = priors_vec.at(16);
  double p_model_3_arg = priors_vec.at(17);
  double wpmean_mean_arg = priors_vec.at(18);
  double wpmean_sd_arg = priors_vec.at(19);
  double dpmean_mean_arg = priors_vec.at(20);
  double dpmean_sd_arg = priors_vec.at(21);
  double wrmean_mean_arg = priors_vec.at(22);
  double wrmean_sd_arg = priors_vec.at(23);
  double wpsd_scale_arg = priors_vec.at(24);
  double tpsd_scale_arg = priors_vec.at(25);
  double dpsd_scale_arg = priors_vec.at(26);
  double wrsd_scale_arg = priors_vec.at(27);
  double sigma_scale_arg = priors_vec.at(28);
  
  print_pos(__FILE__,__LINE__);
  
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
  
  print_pos(__FILE__,__LINE__);
  
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
  
  print_pos(__FILE__,__LINE__);
  
  // Initialize likelihood
  current_parameters.log_likelihood = log_likelihood(viral_data,
                                                     current_data,
                                                     current_parameters,
                                                     settings);
  
  print_pos(__FILE__,__LINE__);
  
  // Initialize output
  std::ofstream scalars_out;
  
  scalars_out.open((std::string(output_dir)+std::string("scalars_out.csv")).c_str());
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << "wp_mean_" << st << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << "wp_sd_" << st << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << "dp_mean_" << st << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << "dp_sd_" << st << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << "tp_sd_" << st << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << "wr_mean_" << st << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << "wr_sd_" << st << ",";
  }
  
  scalars_out << "sigma,log_likelihood\n";
  
  // Output initial values
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << current_parameters.wp_mean.at(st) << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << current_parameters.wp_sd.at(st) << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << current_parameters.dp_mean.at(st) << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << current_parameters.dp_sd.at(st) << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << current_parameters.tp_sd.at(st) << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << current_parameters.wr_mean.at(st) << ",";
  }
  for(int st = 0; st < settings.n_subtypes; st++){
    scalars_out << current_parameters.wr_sd.at(st) << ",";
  }
  
  scalars_out << current_parameters.sigma << "," << current_parameters.log_likelihood << "\n";
  
  // Begin loop
  double norm_draw;
  double unif_draw;
  
  for(int iter = 0; iter < settings.n_iterations; iter++){
    // At each iteration:
    // Update wp_mean for each subtype
    for(int st = 0; st < settings.n_subtypes; st++){
      norm_draw = rnorm_boost(0,1,rng_value);
      unif_draw = runif(0,1,rng_value);
      
      update_wp_mean(st,
                     current_data,
                     current_parameters,
                     settings,
                     viral_data,
                     priors,
                     scaling_factors,
                     norm_draw,
                     unif_draw);
    }
    
    // Update wp_sd for each type
    for(int st = 0; st < settings.n_subtypes; st++){
      norm_draw = rnorm_boost(0,1,rng_value);
      unif_draw = runif(0,1,rng_value);
      
      update_wp_sd(st,
                     current_data,
                     current_parameters,
                     settings,
                     viral_data,
                     priors,
                     scaling_factors,
                     norm_draw,
                     unif_draw);
    }
    
    // Update dp_mean for each type
    for(int st = 0; st < settings.n_subtypes; st++){
      norm_draw = rnorm_boost(0,1,rng_value);
      unif_draw = runif(0,1,rng_value);
      
      update_dp_mean(st,
                     current_data,
                     current_parameters,
                     settings,
                     viral_data,
                     priors,
                     scaling_factors,
                     norm_draw,
                     unif_draw);
    }
    
    // Update dp_sd for each type
    for(int st = 0; st < settings.n_subtypes; st++){
      norm_draw = rnorm_boost(0,1,rng_value);
      unif_draw = runif(0,1,rng_value);
      
      update_dp_sd(st,
                     current_data,
                     current_parameters,
                     settings,
                     viral_data,
                     priors,
                     scaling_factors,
                     norm_draw,
                     unif_draw);
    }
    
    // Update tp_sd for each type
    for(int st = 0; st < settings.n_subtypes; st++){
      norm_draw = rnorm_boost(0,1,rng_value);
      unif_draw = runif(0,1,rng_value);
      
      update_tp_sd(st,
                     current_data,
                     current_parameters,
                     settings,
                     viral_data,
                     priors,
                     scaling_factors,
                     norm_draw,
                     unif_draw);
    }
    
    // Update wr_mean for each type
    for(int st = 0; st < settings.n_subtypes; st++){
      norm_draw = rnorm_boost(0,1,rng_value);
      unif_draw = runif(0,1,rng_value);
      
      update_wr_mean(st,
                     current_data,
                     current_parameters,
                     settings,
                     viral_data,
                     priors,
                     scaling_factors,
                     norm_draw,
                     unif_draw);
    }
    
    // Update wr_sd for each type
    for(int st = 0; st < settings.n_subtypes; st++){
      norm_draw = rnorm_boost(0,1,rng_value);
      unif_draw = runif(0,1,rng_value);
      
      update_wr_sd(st,
                     current_data,
                     current_parameters,
                     settings,
                     viral_data,
                     priors,
                     scaling_factors,
                     norm_draw,
                     unif_draw);
    }
    
    norm_draw = rnorm_boost(0,1,rng_value);
    unif_draw = runif(0,1,rng_value);
    
    // Update sigma
    update_sigma(current_data,
                   current_parameters,
                   settings,
                   viral_data,
                   priors,
                   scaling_factors,
                   norm_draw,
                   unif_draw);
    
    // Update model for each individual
    // Update wp for each individual
    // Update dp for each individual
    // Update tp for each individual
    // Update wr for each individual
    // Output results
    // Output scalars
    // Output initial values
    for(int st = 0; st < settings.n_subtypes; st++){
      scalars_out << current_parameters.wp_mean.at(st) << ",";
    }
    for(int st = 0; st < settings.n_subtypes; st++){
      scalars_out << current_parameters.wp_sd.at(st) << ",";
    }
    for(int st = 0; st < settings.n_subtypes; st++){
      scalars_out << current_parameters.dp_mean.at(st) << ",";
    }
    for(int st = 0; st < settings.n_subtypes; st++){
      scalars_out << current_parameters.dp_sd.at(st) << ",";
    }
    for(int st = 0; st < settings.n_subtypes; st++){
      scalars_out << current_parameters.tp_sd.at(st) << ",";
    }
    for(int st = 0; st < settings.n_subtypes; st++){
      scalars_out << current_parameters.wr_mean.at(st) << ",";
    }
    for(int st = 0; st < settings.n_subtypes; st++){
      scalars_out << current_parameters.wr_sd.at(st) << ",";
    }
    
    scalars_out << current_parameters.sigma << "," << current_parameters.log_likelihood << "\n";
  
    Rcout << "Iteration: " << iter << "\n";
  }
  
}