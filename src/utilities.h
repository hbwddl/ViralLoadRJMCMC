#ifndef UTILITIES_H
#define UTILITIES_H

#include "structs.h"

double mu(double t, double wp, double tp, double dp, double wr);
void print_settings_struct(settings_struct& settings_struct_arg);
void print_double_vec(std::vector<double> double_vec);
void print_int_vec(std::vector<int> int_vec);
void print_viral_data_struct(viral_data_struct& viral_data_struct_arg);
void print_current_data_struct(current_data_struct& current_data_struct_arg);
void print_current_parameters(current_parameters_struct& current_parameters_struct_arg);
void print_priors(priors_struct& priors_struct_arg);
              
#endif