#ifndef LOG_LIKELIHOOD_H
#define LOG_LIKELIHOOD_H

#include <boost/math/distributions/normal.hpp>
#include "structs.h"

double log_likelihood_ti(double y_i, double t_i, double wp_i, double tp_i, double dp_i, double wr_i, double sigma, double sensitivity);
double log_likelihood(viral_data_struct& viral_data,
                      current_data_struct& current_data,
                      current_parameters_struct& current_parameters,
                      settings_struct& settings);

#endif