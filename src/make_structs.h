#ifndef MAKE_STRUCTS_H
#define MAKE_STRUCTS_H

#include "structs.h"

void make_settings_struct(settings_struct& settings_struct_arg,
                          double lod,
                          double sensitivity,
                          double n_iterations,
                          double n_subtypes,
                          double n_subjects,
                          double n_data);

void make_viral_data_struct(viral_data_struct& viral_data_struct_arg,
                            std::vector<int> index, // length n_data
                            std::vector<double> viral_load, 
                            std::vector<double> time, 
                            std::vector<int> subtype, // length n_subjects
                            std::vector<double> t_first_positive,
                            std::vector<double> t_last_positive,
                            std::vector<double> t_first_test,
                            std::vector<double> t_last_test,
                            std::vector<int> n_positive_tests);

void make_current_data_struct(current_data_struct& current_data_struct_arg,
                              std::vector<double> wp_current,
                              std::vector<double> tp_current,
                              std::vector<double> dp_current,
                              std::vector<double> wr_current,
                              std::vector<int> model_current);

void make_current_parameters_struct(current_parameters_struct& current_parameters_struct_arg,
                               double log_likelihood,
                               double sigma,
                               std::vector<double> wp_mean,
                               std::vector<double> wp_sd,
                               std::vector<double> tp_sd,
                               std::vector<double> dp_mean,
                               std::vector<double> dp_sd,
                               std::vector<double> wr_mean,
                               std::vector<double> wr_sd);

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
                        double sigma_scale);

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
                                 double wr_sf);


#endif