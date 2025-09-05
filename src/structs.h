#ifndef STRUCTS_H
#define STRUCTS_H

struct settings_struct{
  double lod; // limit of detection.
  double sensitivity; // false positive rate
  double n_iterations; // number of iterations
  double n_subtypes; // number of subtypes comparing
  double n_subjects; // number of subjects in analysis
  double n_data; // number of datapoints in analysis
};

struct viral_data_struct{
  // Will be length n_data
  std::vector<int> index;
  std::vector<double> viral_load;
  std::vector<double> time;
  
  // Will be length n_subjects
  std::vector<int> subtype;
  std::vector<double> t_first_positive;
  std::vector<double> t_last_positive;
  std::vector<double> t_first_test;
  std::vector<double> t_last_test;
  std::vector<int> n_positive_tests;
  std::vector<double> max_viral_load;
  
};

struct current_data_struct{
  std::vector<double> wp_current;
  std::vector<double> tp_current;
  std::vector<double> dp_current;
  std::vector<double> wr_current;
  std::vector<int> model_current;
};

struct current_parameters_struct{
  double log_likelihood;
  double sigma;
  std::vector<double> wp_mean;
  std::vector<double> wp_sd;
  std::vector<double> tp_sd;
  std::vector<double> dp_mean;
  std::vector<double> dp_sd;
  std::vector<double> wr_mean;
  std::vector<double> wr_sd;
};

struct priors_struct{
  double wp_min;
  double wp_max;
  double wr_min;
  double wr_max;
  double wpmean_max;
  double dpmean_max;
  double wrmean_max;
  double wpsd_max;
  double tpsd_max;
  double dpsd_max;
  double wrsd_max;
  double sigma_max;
  double wpsd_min;
  double tpsd_min;
  double dpsd_min;
  double wrsd_min;
  double sigma_min;
  double p_model_1;
  double p_model_2;
  double p_model_3;
  
  double wpmean_mean;
  double wpmean_sd;
  double dpmean_mean;
  double dpmean_sd;
  double wrmean_mean;
  double wrmean_sd;
  
  double wpsd_scale;
  double tpsd_scale;
  double dpsd_scale;
  double wrsd_scale;
  double sigma_scale;
  
};

struct scaling_factors_struct{
  std::vector<double> wp_mean_sf;
  std::vector<double> wp_sd_sf;
  std::vector<double> tp_sd_sf;
  std::vector<double> dp_mean_sf;
  std::vector<double> dp_sd_sf;
  std::vector<double> wr_mean_sf;
  std::vector<double> wr_sd_sf;
  
  double sigma_sf;
  
  double wp_sf;
  double tp_sf;
  double dp_sf;
  double wr_sf;
  
};

#endif