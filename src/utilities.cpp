#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;

// [[Rcpp::export]]
double mu(double t, double wp, double tp, double dp, double wr){
  if(t > (tp - wp) && t <= tp){
    return((dp/wp)*(t-(tp-wp)));
  } else if(t > tp && t <= (tp + wr)){
    return(dp - (dp/wr)*(t - tp));
  } else{
    return(0.0);
  }
}

void print_settings_struct(settings_struct& settings_struct_arg){
  Rcout << "Settings:" << "\n";
  Rcout << "L.O.D.: " << settings_struct_arg.lod << "\n";
  Rcout << "Sensitivity: " << settings_struct_arg.sensitivity << "\n";
  Rcout << "N Iterations: " << settings_struct_arg.n_iterations << "\n";
  Rcout << "N Subtypes: " << settings_struct_arg.n_subtypes << "\n";
  Rcout << "N Subjects: " << settings_struct_arg.n_subjects << "\n";
  Rcout << "N Data: " << settings_struct_arg.n_data << "\n";
}

void print_double_vec(std::vector<double> double_vec){
  Rcout << "Size: " << double_vec.size() << "\n";
  for(int i = 0; i < double_vec.size(); i++){
    Rcout << double_vec.at(i) << ",";
  }
  Rcout << "\n";
}

void print_int_vec(std::vector<int> int_vec){
  Rcout << "Size: " << int_vec.size() << "\n";
  for(int i = 0; i < int_vec.size(); i++){
    Rcout << int_vec.at(i) << ",";
  }
  Rcout << "\n";
}

void print_viral_data_struct(viral_data_struct& viral_data_struct_arg){
  Rcout << "Viral data:\n";
  Rcout << "Index:\n";
  print_int_vec(viral_data_struct_arg.index);
  
  Rcout << "Viral Load:\n";
  print_double_vec(viral_data_struct_arg.viral_load);
  
  Rcout << "Time:\n";
  print_double_vec(viral_data_struct_arg.time);
  
  Rcout << "Subtype:\n";
  print_int_vec(viral_data_struct_arg.subtype);
  
  Rcout << "First positive t:\n";
  print_double_vec(viral_data_struct_arg.t_first_positive);
  
  Rcout << "Last positive t:\n";
  print_double_vec(viral_data_struct_arg.t_last_positive);
  
  Rcout << "First test t:\n";
  print_double_vec(viral_data_struct_arg.t_first_test);
  
  Rcout << "Last test t:\n";
  print_double_vec(viral_data_struct_arg.t_last_test);
  
  Rcout << "N positive tests:\n";
  print_int_vec(viral_data_struct_arg.n_positive_tests);
}

void print_current_data_struct(current_data_struct& current_data_struct_arg){
  Rcout << "Current data:\n";
  Rcout << "WP Current:\n";
  print_double_vec(current_data_struct_arg.wp_current);
  
  Rcout << "TP Current:\n";
  print_double_vec(current_data_struct_arg.tp_current);
  
  Rcout << "DP Current:\n";
  print_double_vec(current_data_struct_arg.dp_current);
  
  Rcout << "WR Current:\n";
  print_double_vec(current_data_struct_arg.wr_current);
  
  Rcout << "Model current:\n";
  print_int_vec(current_data_struct_arg.model_current);
}

void print_current_parameters(current_parameters_struct& current_parameters_struct_arg){
  Rcout << "Current parameters:";
  Rcout << "Likelihood: " << current_parameters_struct_arg.log_likelihood << "\n";
  Rcout << "Sigma: " << current_parameters_struct_arg.sigma << "\n";
  
  Rcout << "WP Mean:\n";
  print_double_vec(current_parameters_struct_arg.wp_mean);
  
  Rcout << "WP SD:\n";
  print_double_vec(current_parameters_struct_arg.wp_sd);
  
  Rcout << "TP SD:\n";
  print_double_vec(current_parameters_struct_arg.tp_sd);
  
  Rcout << "DP Mean:\n";
  print_double_vec(current_parameters_struct_arg.dp_mean);
  
  Rcout << "DP SD:\n";
  print_double_vec(current_parameters_struct_arg.dp_sd);
  
  Rcout << "WR Mean:\n";
  print_double_vec(current_parameters_struct_arg.wr_mean);
  
  Rcout << "WR SD:\n";
  print_double_vec(current_parameters_struct_arg.wr_sd);
}

void print_priors(priors_struct& priors_struct_arg){
  Rcout << "Priors:\n";
  Rcout << "WP Max: " << priors_struct_arg.wp_max << "\n";
  Rcout << "WR Max: " << priors_struct_arg.wr_max << "\n";
  Rcout << "WP Mean Max: " << priors_struct_arg.wpmean_max << "\n";
  Rcout << "DP Mean Max: " << priors_struct_arg.dpmean_max << "\n";
  Rcout << "WR Mean Max: " << priors_struct_arg.wrmean_max << "\n";
  Rcout << "WP SD Max: " << priors_struct_arg.wpsd_max << "\n";
  Rcout << "TP SD Max: " << priors_struct_arg.tpsd_max << "\n";
  Rcout << "DP SD Max: " << priors_struct_arg.dpsd_max << "\n";
  Rcout << "WR SD Max: " << priors_struct_arg.wrsd_max << "\n";
  Rcout << "Sigma Max: " << priors_struct_arg.sigma_max << "\n";
  Rcout << "WP SD Min: " << priors_struct_arg.wpsd_min << "\n";
  Rcout << "TP SD Min: " << priors_struct_arg.tpsd_min << "\n";
  Rcout << "DP SD Min: " << priors_struct_arg.dpsd_min << "\n";
  Rcout << "WR SD Min: " << priors_struct_arg.wrsd_min << "\n";
  Rcout << "Sigma Min: " << priors_struct_arg.sigma_min << "\n";
  Rcout << "Model P: " << priors_struct_arg.p_model_1 << " " << priors_struct_arg.p_model_2 << " " << priors_struct_arg.p_model_3 << "\n";
  
  Rcout << "WP Mean Mean: " << priors_struct_arg.wpmean_mean << "\n";
  Rcout << "WP Mean SD: " << priors_struct_arg.wpmean_sd << "\n";
  Rcout << "DP Mean Mean: " << priors_struct_arg.dpmean_mean << "\n";
  Rcout << "DP Mean SD: " << priors_struct_arg.dpmean_sd << "\n";
  Rcout << "WR Mean Mean: " << priors_struct_arg.wrmean_mean << "\n";
  Rcout << "WR Mean SD: " << priors_struct_arg.wrmean_sd << "\n";
  
  Rcout << "WP SD Scale: " << priors_struct_arg.wpsd_scale << "\n";
  Rcout << "TP SD Scale: " << priors_struct_arg.tpsd_scale << "\n";
  Rcout << "DP SD Scale: " << priors_struct_arg.dpsd_scale << "\n";
  Rcout << "WR SD Scale: " << priors_struct_arg.wrsd_scale << "\n";
  Rcout << "Sigma Scale: " << priors_struct_arg.sigma_scale << "\n";
}