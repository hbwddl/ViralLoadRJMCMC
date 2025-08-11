#include <Rcpp.h>
#include "make_structs.h"

using namespace Rcpp;

void make_settings_struct(settings_struct& settings_struct_arg,
                          double lod,
                          double sensitivity,
                          double n_iterations,
                          double n_subtypes,
                          double n_subjects,
                          double n_data){
  settings_struct_arg.lod = lod; 
  settings_struct_arg.sensitivity = sensitivity; 
  settings_struct_arg.n_iterations = n_iterations; 
  settings_struct_arg.n_subtypes = n_subtypes; 
  settings_struct_arg.n_subjects = n_subjects; 
  settings_struct_arg.n_data = n_data; 
}

void make_viral_data_struct(viral_data_struct& viral_data_struct_arg,
                            std::vector<int> index, // length n_data
                            std::vector<double> viral_load, 
                            std::vector<double> time, 
                            std::vector<int> subtype, // length n_subjects
                            std::vector<double> t_first_positive,
                            std::vector<double> t_last_positive,
                            std::vector<double> t_first_test,
                            std::vector<double> t_last_test,
                            std::vector<int> n_positive_tests
                            ){
  viral_data_struct_arg.index = index;
  viral_data_struct_arg.viral_load = viral_load;
  viral_data_struct_arg.time = time;
  viral_data_struct_arg.subtype = subtype;
  viral_data_struct_arg.t_first_positive = t_first_positive;
  viral_data_struct_arg.t_last_positive = t_last_positive;
  viral_data_struct_arg.t_first_test = t_first_test;
  viral_data_struct_arg.t_last_test = t_last_test;
  viral_data_struct_arg.n_positive_tests = n_positive_tests;
}

void make_current_data_struct(current_data_struct& current_data_struct_arg,
                              std::vector<double> wp_current,
                              std::vector<double> tp_current,
                              std::vector<double> dp_current,
                              std::vector<double> wr_current,
                              std::vector<int> model_current){
  current_data_struct_arg.wp_current = wp_current;
  current_data_struct_arg.tp_current = tp_current;
  current_data_struct_arg.dp_current = dp_current;
  current_data_struct_arg.wr_current = wr_current;
  current_data_struct_arg.model_current = model_current;
}

void make_current_parameters_struct(current_parameters_struct& current_parameters_struct_arg,
                                    double log_likelihood,
                                    double sigma,
                                    std::vector<double> wp_mean,
                                    std::vector<double> wp_sd,
                                    std::vector<double> tp_sd,
                                    std::vector<double> dp_mean,
                                    std::vector<double> dp_sd,
                                    std::vector<double> wr_mean,
                                    std::vector<double> wr_sd){
  
  current_parameters_struct_arg.log_likelihood = log_likelihood;
  current_parameters_struct_arg.sigma = sigma;
  current_parameters_struct_arg.wp_mean = wp_mean;
  current_parameters_struct_arg.wp_sd = wp_sd;
  current_parameters_struct_arg.tp_sd = tp_sd;
  current_parameters_struct_arg.dp_mean = dp_mean;
  current_parameters_struct_arg.dp_sd = dp_sd;
  current_parameters_struct_arg.wr_mean = wr_mean;
  current_parameters_struct_arg.wr_sd = wr_sd;
  
}

void make_priors_struct(priors_struct& priors_struct_arg,
                        double wp_max,
                        double wr_max,
                        double wpmean_max,
                        double dpmean_max,
                        double wrmean_max,
                        double wpsd_max,
                        double tpsd_max,
                        double dpsd_max,
                        double wrsd_max,
                        double sigma_max,
                        double wpsd_min,
                        double tpsd_min,
                        double dpsd_min,
                        double wrsd_min,
                        double sigma_min,
                        double p_model_1,
                        double p_model_2,
                        double p_model_3,
                        double wpmean_mean,
                        double wpmean_sd,
                        double dpmean_mean,
                        double dpmean_sd,
                        double wrmean_mean,
                        double wrmean_sd,
                        double wpsd_scale,
                        double tpsd_scale,
                        double dpsd_scale,
                        double wrsd_scale,
                        double sigma_scale){
  priors_struct_arg.wp_max = wp_max;
  priors_struct_arg.wr_max = wr_max;
  priors_struct_arg.wpmean_max = wpmean_max;
  priors_struct_arg.dpmean_max = dpmean_max;
  priors_struct_arg.wrmean_max = wrmean_max;
  priors_struct_arg.wpsd_max = wpsd_max;
  priors_struct_arg.tpsd_max = tpsd_max;
  priors_struct_arg.dpsd_max = dpsd_max;
  priors_struct_arg.wrsd_max = wrsd_max;
  priors_struct_arg.sigma_max = sigma_max;
  priors_struct_arg.wpsd_min = wpsd_min;
  priors_struct_arg.tpsd_min = tpsd_min;
  priors_struct_arg.dpsd_min = dpsd_min;
  priors_struct_arg.wrsd_min = wrsd_min;
  priors_struct_arg.sigma_min = sigma_min;
  priors_struct_arg.p_model_1 = p_model_1;
  priors_struct_arg.p_model_2 = p_model_2;
  priors_struct_arg.p_model_3 = p_model_3;
  priors_struct_arg.wpmean_mean = wpmean_mean;
  priors_struct_arg.wpmean_sd = wpmean_sd;
  priors_struct_arg.dpmean_mean = dpmean_mean;
  priors_struct_arg.dpmean_sd = dpmean_sd;
  priors_struct_arg.wrmean_mean = wrmean_mean;
  priors_struct_arg.wrmean_sd = wrmean_sd;
  priors_struct_arg.wpsd_scale = wpsd_scale;
  priors_struct_arg.tpsd_scale = tpsd_scale;
  priors_struct_arg.dpsd_scale = dpsd_scale;
  priors_struct_arg.wrsd_scale = wrsd_scale;
  priors_struct_arg.sigma_scale = sigma_scale;
}

void make_scaling_factors_struct(scaling_factors_struct& scaling_factors_struct_arg,
                                 std::vector<double> wp_mean_sf,
                                 std::vector<double> wp_sd_sf,
                                 std::vector<double> tp_sd_sf,
                                 std::vector<double> dp_mean_sf,
                                 std::vector<double> dp_sd_sf,
                                 std::vector<double> wr_mean_sf,
                                 std::vector<double> wr_sd_sf,
                                 double sigma_sf,
                                 double wp_sf,
                                 double tp_sf,
                                 double dp_sf,
                                 double wr_sf){
  scaling_factors_struct_arg.wp_mean_sf = wp_mean_sf;
  scaling_factors_struct_arg.wp_sd_sf = wp_sd_sf;
  scaling_factors_struct_arg.tp_sd_sf = tp_sd_sf;
  scaling_factors_struct_arg.dp_mean_sf = dp_mean_sf;
  scaling_factors_struct_arg.dp_sd_sf = dp_sd_sf;
  scaling_factors_struct_arg.wr_mean_sf = wr_mean_sf;
  scaling_factors_struct_arg.wr_sd_sf = wr_sd_sf;
  scaling_factors_struct_arg.sigma_sf = sigma_sf;
  scaling_factors_struct_arg.wp_sf = wp_sf;
  scaling_factors_struct_arg.tp_sf = tp_sf;
  scaling_factors_struct_arg.dp_sf = dp_sf;
  scaling_factors_struct_arg.wr_sf = wr_sf;
}