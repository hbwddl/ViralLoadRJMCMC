#ifndef UPDATE_MODEL_H
#define UPDATE_MODEL_H

#include "structs.h"

void update_model_i(int index_update,
                    current_data_struct& current_data_arg,
                    current_parameters_struct& current_parameters_arg,
                    settings_struct& settings_arg,
                    viral_data_struct& viral_data_arg,
                    priors_struct& priors_arg,
                    scaling_factors_struct& scaling_factors_arg,
                    double norm_0_1_draw,
                    double unif_0_1_draw_model,
                    double unif_0_1_draw_rjmcmc_u,
                    double unif_0_1_draw_rjmcmc_v,
                    double unif_0_1_draw_acp);

#endif