#ifndef UPDATE_INDIVIDUAL_DATA_H
#define UPDATE_INDIVIDUAL_DATA_H

#include "structs.h"

void update_wp_i(int index_update,
                 current_data_struct& current_data_arg,
                 current_parameters_struct& current_parameters_arg,
                 settings_struct& settings_arg,
                 viral_data_struct& viral_data_arg,
                 priors_struct& priors_arg,
                 scaling_factors_struct& scaling_factors_arg,
                 double norm_0_1_draw,
                 double unif_0_1_draw);

void update_tp_i(int index_update,
                 current_data_struct& current_data_arg,
                 current_parameters_struct& current_parameters_arg,
                 settings_struct& settings_arg,
                 viral_data_struct& viral_data_arg,
                 priors_struct& priors_arg,
                 scaling_factors_struct& scaling_factors_arg,
                 double norm_0_1_draw,
                 double unif_0_1_draw);

void update_dp_i(int index_update,
                 current_data_struct& current_data_arg,
                 current_parameters_struct& current_parameters_arg,
                 settings_struct& settings_arg,
                 viral_data_struct& viral_data_arg,
                 priors_struct& priors_arg,
                 scaling_factors_struct& scaling_factors_arg,
                 double norm_0_1_draw,
                 double unif_0_1_draw);

void update_wr_i(int index_update,
                 current_data_struct& current_data_arg,
                 current_parameters_struct& current_parameters_arg,
                 settings_struct& settings_arg,
                 viral_data_struct& viral_data_arg,
                 priors_struct& priors_arg,
                 scaling_factors_struct& scaling_factors_arg,
                 double norm_0_1_draw,
                 double unif_0_1_draw);

#endif