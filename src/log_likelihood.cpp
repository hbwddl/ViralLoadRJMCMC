#include <Rcpp.h>
#include "structs.h"
#include "utilities.h"
#include "distributions.h"
#include "log_likelihood.h"
#include "debug.h"

using namespace Rcpp;

#define debug_lh 0

// [[Rcpp::export]]
double log_likelihood_ti(double y_i, double t_i, double wp_i, double tp_i, double dp_i, double wr_i, double sigma, double sensitivity){ // Lambda = test sensitivity
  double mu_i = mu(t_i, wp_i, tp_i, dp_i, wr_i);
  
  double dev = y_i - mu_i;
  
  return(log(sensitivity*pdf_normal(dev, 0.0, sigma) + (1-sensitivity)*pdf_exponential(y_i,1/log(10))));
  
  if(mu_i > 0 & y_i == 0){ // False negative?
    return(log(1-cdf_normal(mu_i, 0.0, sigma)));
  } else{
    return(log(sensitivity*pdf_normal(dev, 0.0, sigma) + (1-sensitivity)*pdf_exponential(y_i,1/log(10))));
  }
  
}

double log_likelihood(viral_data_struct& viral_data,
                      current_data_struct& current_data,
                      current_parameters_struct& current_parameters,
                      settings_struct& settings){
  long double log_lh_total = 0.0;
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  // Rcout << viral_data.viral_load.size() << "\n";
  
  
  for(int i = 0; i < viral_data.viral_load.size(); i++){
    log_lh_total += log_likelihood_ti(viral_data.viral_load.at(i), 
                                      viral_data.time.at(i), 
                                      current_data.wp_current.at(viral_data.index.at(i)), 
                                      current_data.tp_current.at(viral_data.index.at(i)), 
                                      current_data.dp_current.at(viral_data.index.at(i)), 
                                      current_data.wr_current.at(viral_data.index.at(i)), 
                                      current_parameters.sigma, 
                                      settings.sensitivity);
  }
  
  // Rcout << __LINE__ << " " << log_lh_total << "\n";
  
  // Rcout << log_lh_total << "\n";
  print_pos(__FILE__,__LINE__,debug_lh);
  
  int model_i;
  int subtype_i;
  double wp_i;
  double wpmean_subtype_i;
  double wpsd_subtype_i;
  double tp_i;
  double tpsd_subtype_i;
  double dp_i;
  double dpmean_subtype_i;
  double dpsd_subtype_i;
  double wr_i;
  double wrmean_subtype_i;
  double wrsd_subtype_i;
  
  double wp_lh_i;
  double tp_lh_i;
  double dp_lh_i;
  double wr_lh_i;
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  for(int i = 0; i < current_data.wp_current.size(); i++){
      model_i = -1;
      subtype_i = -1;
      wp_i = -1;
      wpmean_subtype_i = -1;
      wpsd_subtype_i = -1;
      
      model_i = current_data.model_current.at(i);
    
      if(model_i == 1 || model_i == 2){
        subtype_i = viral_data.subtype.at(i);
        wp_i = current_data.wp_current.at(i);
        wpmean_subtype_i = current_parameters.wp_mean.at(subtype_i);
        wpsd_subtype_i = current_parameters.wp_sd.at(subtype_i);
        wp_lh_i = 0;
        wp_lh_i = log(pdf_normal(wp_i,
                                 wpmean_subtype_i,
                                 wpsd_subtype_i));
        log_lh_total += wp_lh_i;
        
        // Rcout << i << "," << inf_type[i] << "," << model[i] << "," << wp_lh_i << "\n";
        
      }
    
  }
  
  // Rcout << log_lh_total << "\n";
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  for(int i = 0; i < current_data.tp_current.size(); i++){
    subtype_i = -1;
    tp_i = -1;
    tpsd_subtype_i = -1;
    
    subtype_i = viral_data.subtype.at(i);
    tp_i =  current_data.tp_current.at(i);
    tpsd_subtype_i = current_parameters.tp_sd.at(subtype_i);
    
    log_lh_total += log(pdf_normal(tp_i,
                                   0,
                                   tpsd_subtype_i));
  }
  
  // Rcout << log_lh_total << "\n";
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  for(int i = 0; i < current_data.dp_current.size(); i++){
    subtype_i = -1;
    dp_i = -1;
    dpmean_subtype_i = -1;
    dpsd_subtype_i = -1;
    
    subtype_i = viral_data.subtype.at(i);
    dp_i =  current_data.dp_current.at(i);
    dpmean_subtype_i = current_parameters.dp_mean.at(subtype_i);
    dpsd_subtype_i = current_parameters.dp_sd.at(subtype_i);
    
    log_lh_total += log(pdf_normal(dp_i,
                                   dpmean_subtype_i,
                                   dpsd_subtype_i));
  }
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  for(int i = 0; i < current_data.wr_current.size(); i++){
    model_i = -1;
    subtype_i = -1;
    wr_i = -1;
    wrmean_subtype_i = -1;
    wrsd_subtype_i = -1;
    
    model_i = current_data.model_current.at(i);
    
    if(model_i == 2 || model_i == 3){
      subtype_i = viral_data.subtype.at(i);
      wr_i = current_data.wr_current.at(i);
      wrmean_subtype_i = current_parameters.wr_mean.at(subtype_i);
      wrsd_subtype_i = current_parameters.wr_sd.at(subtype_i);
      
      log_lh_total += log(pdf_normal(wr_i,
                                 wrmean_subtype_i,
                                 wrsd_subtype_i));
    }
  }
  
  print_pos(__FILE__,__LINE__,debug_lh);
  // Rcout << log_lh_total << "\n";
  
  return(log_lh_total);
}


void check_log_likelihood(double current_likelihood_arg,
                          viral_data_struct& viral_data_arg,
                          current_data_struct& current_data_arg,
                          current_parameters_struct& current_parameters_arg,
                          settings_struct& settings_arg,
                          double tolerance){
  double check_likelihood = log_likelihood(viral_data_arg,
                                           current_data_arg,
                                           current_parameters_arg,
                                           settings_arg);
  
  if(std::abs(check_likelihood - current_likelihood_arg) > tolerance){
    Rcout << "ERR LIKELIHOOD CHECK Current " << current_likelihood_arg << " Calculated " << check_likelihood << "\n";
  }
}

void check_data(viral_data_struct& viral_data_arg,
                current_data_struct& current_data_arg,
                current_parameters_struct& current_parameters_arg,
                settings_struct& settings_arg,
                priors_struct& priors_arg){
  
  for(int subj = 0; subj < settings_arg.n_subjects; subj++){
    /* Check wp, wr, tp, dp within bounds */
    if(current_data_arg.wp_current.at(subj) < priors_arg.wp_min || 
        current_data_arg.wp_current.at(subj) > priors_arg.wp_max){
      Rcout << "ERR WP SUBJ " << subj << " WP: " << current_data_arg.wp_current.at(subj) << "\n";
    }
    if(current_data_arg.dp_current.at(subj) < (viral_data_arg.max_viral_load.at(subj)/2) || 
       current_data_arg.dp_current.at(subj) > settings_arg.lod){
      Rcout << "ERR DP SUBJ " << subj << " DP: " << current_data_arg.dp_current.at(subj) << "\n";
    }
    
    if(current_data_arg.model_current.at(subj) == 2 &&
       ((current_data_arg.tp_current.at(subj) < viral_data_arg.t_first_test.at(subj) || 
       current_data_arg.tp_current.at(subj) > viral_data_arg.t_last_test.at(subj)) ||
       (current_data_arg.tp_current.at(subj) < -priors_arg.wp_max ||
       current_data_arg.tp_current.at(subj) > priors_arg.wr_max))){
      Rcout << "ERR TP RANGE MODEL 2 SUBJ " << subj << " TP: " << current_data_arg.tp_current.at(subj) << "\n";
    }
    
    if(current_data_arg.model_current.at(subj) == 1 &&
       ((current_data_arg.tp_current.at(subj) < viral_data_arg.t_last_test.at(subj)) ||
       (current_data_arg.tp_current.at(subj) < -priors_arg.wp_max ||
       current_data_arg.tp_current.at(subj) > priors_arg.wr_max))){
      Rcout << "ERR TP RANGE MODEL 1 SUBJ " << subj << " TP: " << current_data_arg.tp_current.at(subj) << "\n";
    }
    
    if(current_data_arg.model_current.at(subj) == 3 &&
       ((current_data_arg.tp_current.at(subj) > viral_data_arg.t_first_test.at(subj)) ||
       (current_data_arg.tp_current.at(subj) < -priors_arg.wp_max ||
       current_data_arg.tp_current.at(subj) > priors_arg.wr_max))){
      Rcout << "ERR TP RANGE MODEL 3 SUBJ " << subj << " TP: " << current_data_arg.tp_current.at(subj) << "\n";
    }
    
    if(current_data_arg.wr_current.at(subj) < priors_arg.wr_min || 
       current_data_arg.wr_current.at(subj) > priors_arg.wr_max){
      Rcout << "ERR WR SUBJ " << subj << " WR: " << current_data_arg.wr_current.at(subj) << "\n";
    }
  }
}

int check_data_err(viral_data_struct& viral_data_arg,
                current_data_struct& current_data_arg,
                current_parameters_struct& current_parameters_arg,
                settings_struct& settings_arg,
                priors_struct& priors_arg){
  
  int data_err = 0;
  
  for(int subj = 0; subj < settings_arg.n_subjects; subj++){
    /* Check wp, wr, tp, dp within bounds */
    if(current_data_arg.wp_current.at(subj) < priors_arg.wp_min || 
    current_data_arg.wp_current.at(subj) > priors_arg.wp_max){
      data_err = 1;
    }
    if(current_data_arg.dp_current.at(subj) < (viral_data_arg.max_viral_load.at(subj)/2) || 
       current_data_arg.dp_current.at(subj) > settings_arg.lod){
      data_err = 1;
    }
    
    if(current_data_arg.model_current.at(subj) == 2 &&
       ((current_data_arg.tp_current.at(subj) < viral_data_arg.t_first_test.at(subj) || 
       current_data_arg.tp_current.at(subj) > viral_data_arg.t_last_test.at(subj)) ||
       (current_data_arg.tp_current.at(subj) < -priors_arg.wp_max ||
       current_data_arg.tp_current.at(subj) > priors_arg.wr_max))){
      data_err = 1;
    }
    
    if(current_data_arg.model_current.at(subj) == 1 &&
       ((current_data_arg.tp_current.at(subj) < viral_data_arg.t_last_test.at(subj)) ||
       (current_data_arg.tp_current.at(subj) < -priors_arg.wp_max ||
       current_data_arg.tp_current.at(subj) > priors_arg.wr_max))){
      data_err = 1;
    }
    
    if(current_data_arg.model_current.at(subj) == 3 &&
       ((current_data_arg.tp_current.at(subj) > viral_data_arg.t_first_test.at(subj)) ||
       (current_data_arg.tp_current.at(subj) < -priors_arg.wp_max ||
       current_data_arg.tp_current.at(subj) > priors_arg.wr_max))){
      data_err = 1;
    }
    
    if(current_data_arg.wr_current.at(subj) < priors_arg.wr_min || 
       current_data_arg.wr_current.at(subj) > priors_arg.wr_max){
      data_err = 1;
    }
  }
  
  return(data_err);
}