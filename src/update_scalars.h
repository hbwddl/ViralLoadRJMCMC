#ifndef UPDATE_SCALARS_H
#define UPDATE_SCALARS_H

#include "structs.h"

void update_wp_mean(int which_subtype_update,
                    current_data_struct& current_data_arg,
                    current_parameters_struct& current_parameters_arg,
                    settings_struct& settings_arg,
                    viral_data_struct& viral_data_arg,
                    priors_struct& priors_arg,
                    scaling_factors_struct& scaling_factors_arg,
                    double norm_0_1_draw,
                    double unif_0_1_draw);

void update_wp_sd(int which_subtype_update,
                  current_data_struct& current_data_arg,
                  current_parameters_struct& current_parameters_arg,
                  settings_struct& settings_arg,
                  viral_data_struct& viral_data_arg,
                  priors_struct& priors_arg,
                  scaling_factors_struct& scaling_factors_arg,
                  double norm_0_1_draw,
                  double unif_0_1_draw);

void update_tp_sd(int which_subtype_update,
                  current_data_struct& current_data_arg,
                  current_parameters_struct& current_parameters_arg,
                  settings_struct& settings_arg,
                  viral_data_struct& viral_data_arg,
                  priors_struct& priors_arg,
                  scaling_factors_struct& scaling_factors_arg,
                  double norm_0_1_draw,
                  double unif_0_1_draw);

void update_dp_mean(int which_subtype_update,
                    current_data_struct& current_data_arg,
                    current_parameters_struct& current_parameters_arg,
                    settings_struct& settings_arg,
                    viral_data_struct& viral_data_arg,
                    priors_struct& priors_arg,
                    scaling_factors_struct& scaling_factors_arg,
                    double norm_0_1_draw,
                    double unif_0_1_draw);

void update_dp_sd(int which_subtype_update,
                  current_data_struct& current_data_arg,
                  current_parameters_struct& current_parameters_arg,
                  settings_struct& settings_arg,
                  viral_data_struct& viral_data_arg,
                  priors_struct& priors_arg,
                  scaling_factors_struct& scaling_factors_arg,
                  double norm_0_1_draw,
                  double unif_0_1_draw);

void update_wr_mean(int which_subtype_update,
                    current_data_struct& current_data_arg,
                    current_parameters_struct& current_parameters_arg,
                    settings_struct& settings_arg,
                    viral_data_struct& viral_data_arg,
                    priors_struct& priors_arg,
                    scaling_factors_struct& scaling_factors_arg,
                    double norm_0_1_draw,
                    double unif_0_1_draw);

void update_wr_sd(int which_subtype_update,
                  current_data_struct& current_data_arg,
                  current_parameters_struct& current_parameters_arg,
                  settings_struct& settings_arg,
                  viral_data_struct& viral_data_arg,
                  priors_struct& priors_arg,
                  scaling_factors_struct& scaling_factors_arg,
                  double norm_0_1_draw,
                  double unif_0_1_draw);

void update_sigma(current_data_struct& current_data_arg,
                  current_parameters_struct& current_parameters_arg,
                  settings_struct& settings_arg,
                  viral_data_struct& viral_data_arg,
                  priors_struct& priors_arg,
                  scaling_factors_struct& scaling_factors_arg,
                  double norm_0_1_draw,
                  double unif_0_1_draw);

#endif