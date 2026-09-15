#include <Rcpp.h>
#include "structs.h"
#include "utilities.h"
#include "distributions.h"
#include "log_likelihood.h"
#include "debug.h"
#include "make_structs.h"

using namespace Rcpp;

#define debug_lh 0
#define shared_sd_lh 1

// [[Rcpp::export]]
double log_likelihood_ti(double y_i, double t_i, double wp_i, double tp_i, double dp_i, double wr_i, double sigma, double sensitivity){ // Lambda = test sensitivity
  double mu_i = mu(t_i, wp_i, tp_i, dp_i, wr_i);
  
  double dev = y_i - mu_i;
  
  // return(log(pdf_normal(dev, 0.0, sigma)));
  
  if(mu_i == 0 & y_i == 0){
    return(0.0);
  }
  
  // return(log(sensitivity*pdf_normal(dev, 0.0, sigma) +
  //           (1-sensitivity)*pdf_exponential(y_i,1/log(10))));

  if(mu_i > 0 & y_i == 0){ // False negative?
    if(std::isinf(log(1-cdf_normal(mu_i, 0.0, sigma)))){
      return(0.0);
    } else{
      return(log(1-cdf_normal(mu_i, 0.0, sigma)));
    }
  } else{
    return(log(sensitivity*pdf_normal(dev, 0.0, sigma) + (1-sensitivity)*pdf_exponential(y_i,1/log(10))));
  }
  
}

double log_likelihood_subject(int subj,
                              viral_data_struct& viral_data,
                              current_data_struct& current_data,
                              current_parameters_struct& current_parameters,
                              settings_struct& settings){
  long double log_lh_total = 0.0;
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  // Rcout << viral_data.viral_load.size() << "\n";
  
  // SHARED SD
  int n_subtype = current_parameters.wp_sd.size();
  
  if(shared_sd_lh == 1){
    if(n_subtype > 1){
      for(int i = 1; i < n_subtype; i++){
        current_parameters.wp_sd.at(i) = current_parameters.wp_sd.at(0);
        current_parameters.dp_sd.at(i) = current_parameters.dp_sd.at(0);
        current_parameters.wr_sd.at(i) = current_parameters.wr_sd.at(0);
      }
    }
  }

  for(int i = 0; i < viral_data.viral_load.size(); i++){
    if(viral_data.index.at(i) == subj){
      log_lh_total += log_likelihood_ti(viral_data.viral_load.at(i), 
                                        viral_data.time.at(i), 
                                        current_data.wp_current.at(viral_data.index.at(i)), 
                                        current_data.tp_current.at(viral_data.index.at(i)), 
                                        current_data.dp_current.at(viral_data.index.at(i)), 
                                        current_data.wr_current.at(viral_data.index.at(i)), 
                                        current_parameters.sigma, 
                                        settings.sensitivity);
    }
    if(debug_lh == 1){
      Rcout << "Line " << i+1 << " subject " << viral_data.index.at(i) << " obs " << viral_data.viral_load.at(i) << 
        " mu " << mu(viral_data.time.at(i), 
                 current_data.wp_current.at(viral_data.index.at(i)), 
                 current_data.tp_current.at(viral_data.index.at(i)), 
                 current_data.dp_current.at(viral_data.index.at(i)), 
                 current_data.wr_current.at(viral_data.index.at(i))) << 
                   " lh_i " << log_likelihood_ti(viral_data.viral_load.at(i), 
                                           viral_data.time.at(i), 
                                           current_data.wp_current.at(viral_data.index.at(i)), 
                                           current_data.tp_current.at(viral_data.index.at(i)), 
                                           current_data.dp_current.at(viral_data.index.at(i)), 
                                           current_data.wr_current.at(viral_data.index.at(i)), 
                                           current_parameters.sigma, 
                                           settings.sensitivity) << "\n";
    }
    
  }
  
  if(debug_lh == 1){
    Rcout << __LINE__ << " " << log_lh_total << "\n";
  }
  
  
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
  
  int i = subj;
  
  model_i = -1;
  subtype_i = -1;
  wp_i = -1;
  wpmean_subtype_i = -1;
  wpsd_subtype_i = -1;
  
  model_i = current_data.model_current.at(subj);
  
  if(model_i == 1 || model_i == 2){
    subtype_i = viral_data.subtype.at(subj);
    wp_i = current_data.wp_current.at(subj);
    wpmean_subtype_i = current_parameters.wp_mean.at(subtype_i);
    wpsd_subtype_i = current_parameters.wp_sd.at(subtype_i);
    wp_lh_i = 0;
    wp_lh_i = log(pdf_normal(wp_i,
                             wpmean_subtype_i,
                             wpsd_subtype_i));
    log_lh_total += wp_lh_i;
    
    if(debug_lh == 1){
      Rcout << subj << "," << subtype_i << "," << model_i << ", wp_lh_i " << wp_lh_i << "\n";
    }
    
  }
  
  if(debug_lh == 1){
    Rcout << log_lh_total << "\n";
  }
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  subtype_i = -1;
  tp_i = -1;
  tpsd_subtype_i = -1;
  
  subtype_i = viral_data.subtype.at(subj);
  tp_i =  current_data.tp_current.at(subj);
  tpsd_subtype_i = current_parameters.tp_sd.at(subtype_i);
  
  tp_lh_i = 0;
  tp_lh_i = log(pdf_normal(tp_i,
                           0,
                           tpsd_subtype_i));
  log_lh_total += tp_lh_i;
  
  if(debug_lh == 1){
    Rcout << subj << "," << subtype_i << "," << model_i << ", tp_i " << tp_i << ", tp_lh_i " << tp_lh_i << "\n";
  }
  
  if(debug_lh == 1){
    Rcout << log_lh_total << "\n";
  }
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  subtype_i = -1;
  dp_i = -1;
  dpmean_subtype_i = -1;
  dpsd_subtype_i = -1;
  
  subtype_i = viral_data.subtype.at(subj);
  dp_i =  current_data.dp_current.at(subj);
  dpmean_subtype_i = current_parameters.dp_mean.at(subtype_i);
  dpsd_subtype_i = current_parameters.dp_sd.at(subtype_i);
  
  dp_lh_i = 0;
  dp_lh_i = log(pdf_normal(dp_i,
                           dpmean_subtype_i,
                           dpsd_subtype_i));
  log_lh_total += dp_lh_i;
  
  if(debug_lh == 1){
    Rcout << subj << "," << subtype_i << "," << model_i << ", dp_lh_i " << dp_lh_i << "\n";
  }
  
  
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  model_i = -1;
  subtype_i = -1;
  wr_i = -1;
  wrmean_subtype_i = -1;
  wrsd_subtype_i = -1;
  
  model_i = current_data.model_current.at(subj);
  
  if(model_i == 2 || model_i == 3){
    subtype_i = viral_data.subtype.at(subj);
    wr_i = current_data.wr_current.at(subj);
    wrmean_subtype_i = current_parameters.wr_mean.at(subtype_i);
    wrsd_subtype_i = current_parameters.wr_sd.at(subtype_i);
    
    wr_lh_i = log(pdf_normal(wr_i,
                             wrmean_subtype_i,
                             wrsd_subtype_i));
    
    log_lh_total += wr_lh_i;
    
    if(debug_lh == 1){
      Rcout << i << "," << subtype_i << "," << model_i << ", wr_lh_i " << wr_lh_i << "\n";
    }
    
  }
  
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  if(debug_lh == 1){
    Rcout << log_lh_total << "\n";
  }
  
  return(log_lh_total);
}

double log_likelihood(viral_data_struct& viral_data,
                      current_data_struct& current_data,
                      current_parameters_struct& current_parameters,
                      settings_struct& settings){
  long double log_lh_total = 0.0;
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  // Rcout << viral_data.viral_load.size() << "\n";
  // SHARED SD
  int n_subtype = current_parameters.wp_sd.size();
  
  if(shared_sd_lh == 1){
    if(n_subtype > 1){
      for(int i = 1; i < n_subtype; i++){
        current_parameters.wp_sd.at(i) = current_parameters.wp_sd.at(0);
        current_parameters.tp_sd.at(i) = current_parameters.tp_sd.at(0);
        current_parameters.dp_sd.at(i) = current_parameters.dp_sd.at(0);
        current_parameters.wr_sd.at(i) = current_parameters.wr_sd.at(0);
      }
    }
  }

  for(int i = 0; i < viral_data.viral_load.size(); i++){
    log_lh_total += log_likelihood_ti(viral_data.viral_load.at(i), 
                                      viral_data.time.at(i), 
                                      current_data.wp_current.at(viral_data.index.at(i)), 
                                      current_data.tp_current.at(viral_data.index.at(i)), 
                                      current_data.dp_current.at(viral_data.index.at(i)), 
                                      current_data.wr_current.at(viral_data.index.at(i)), 
                                      current_parameters.sigma, 
                                      settings.sensitivity);
    
    if(debug_lh == 1){
      Rcout << "Line " << i+1 << " subject " << viral_data.index.at(i) << " obs " << viral_data.viral_load.at(i) << 
        " mu " << mu(viral_data.time.at(i), 
                 current_data.wp_current.at(viral_data.index.at(i)), 
                 current_data.tp_current.at(viral_data.index.at(i)), 
                 current_data.dp_current.at(viral_data.index.at(i)), 
                 current_data.wr_current.at(viral_data.index.at(i))) << 
                   " lh_i " << log_likelihood_ti(viral_data.viral_load.at(i), 
                                           viral_data.time.at(i), 
                                           current_data.wp_current.at(viral_data.index.at(i)), 
                                           current_data.tp_current.at(viral_data.index.at(i)), 
                                           current_data.dp_current.at(viral_data.index.at(i)), 
                                           current_data.wr_current.at(viral_data.index.at(i)), 
                                           current_parameters.sigma, 
                                           settings.sensitivity) << "\n";
    }
    
  }
  
  if(debug_lh == 1){
    Rcout << __LINE__ << " " << log_lh_total << "\n";
  }
  
  
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
      
      if(debug_lh == 1){
        Rcout << i << "," << subtype_i << "," << model_i << ", wp_lh_i " << wp_lh_i << "\n";
      }
      
    }
    
  }
  
  if(debug_lh == 1){
    Rcout << log_lh_total << "\n";
  }
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  for(int i = 0; i < current_data.tp_current.size(); i++){
    subtype_i = -1;
    tp_i = -1;
    tpsd_subtype_i = -1;
    model_i = -1;
    
    model_i = current_data.model_current.at(i);
    subtype_i = viral_data.subtype.at(i);
    tp_i =  current_data.tp_current.at(i);
    tpsd_subtype_i = current_parameters.tp_sd.at(subtype_i);
    
    tp_lh_i = 0;
    tp_lh_i = log(pdf_normal(tp_i,
                             0,
                             tpsd_subtype_i));
    log_lh_total += tp_lh_i;
    
    if(debug_lh == 1){
      Rcout << i << "," << subtype_i << "," << model_i << ", tp_i " << tp_i << ", tp_lh_i " << tp_lh_i << "\n";
    }
    
  }
  
  if(debug_lh == 1){
    Rcout << log_lh_total << "\n";
  }
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  for(int i = 0; i < current_data.dp_current.size(); i++){
    subtype_i = -1;
    dp_i = -1;
    dpmean_subtype_i = -1;
    dpsd_subtype_i = -1;
    model_i = -1;
    
    model_i = current_data.model_current.at(i);
    subtype_i = viral_data.subtype.at(i);
    dp_i =  current_data.dp_current.at(i);
    dpmean_subtype_i = current_parameters.dp_mean.at(subtype_i);
    dpsd_subtype_i = current_parameters.dp_sd.at(subtype_i);
    
    dp_lh_i = 0;
    dp_lh_i = log(pdf_normal(dp_i,
                             dpmean_subtype_i,
                             dpsd_subtype_i));
    log_lh_total += dp_lh_i;
    
    if(debug_lh == 1){
      Rcout << i << "," << subtype_i << "," << model_i << " dp_mean " << dpmean_subtype_i << " dp_sd " << dpsd_subtype_i << ", dp " << dp_i << ", dp_lh_i " << dp_lh_i << "\n";
    }
    
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
      
      wr_lh_i = log(pdf_normal(wr_i,
                               wrmean_subtype_i,
                               wrsd_subtype_i));
      
      log_lh_total += wr_lh_i;
      
      if(debug_lh == 1){
        Rcout << i << "," << subtype_i << "," << model_i << ", wr_lh_i " << wr_lh_i << "\n";
      }
      
    }
  }
  
  print_pos(__FILE__,__LINE__,debug_lh);
  
  if(debug_lh == 1){
    Rcout << log_lh_total << "\n";
  }
  
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
    if(current_data_arg.dp_current.at(subj) < 0 || 
       current_data_arg.dp_current.at(subj) > settings_arg.lod){
      Rcout << "ERR DP SUBJ " << subj << " DP: " << current_data_arg.dp_current.at(subj) << "\n";
    }
    
    if(current_data_arg.model_current.at(subj) == 2 &&
       ((current_data_arg.tp_current.at(subj) < viral_data_arg.t_first_test.at(subj) || 
       current_data_arg.tp_current.at(subj) > viral_data_arg.t_last_test.at(subj)))){
      Rcout << "ERR TP RANGE MODEL 2 SUBJ " << subj << " TP: " << current_data_arg.tp_current.at(subj) << " " << 
        viral_data_arg.t_first_test.at(subj) << " " << 
          viral_data_arg.t_last_test.at(subj) << " " <<
            priors_arg.wp_max << " " <<
              priors_arg.wr_max << "\n";
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
      Rcout << "ERR TP RANGE MODEL 3 SUBJ " << subj << " TP: " << current_data_arg.tp_current.at(subj) << " first positive " << viral_data_arg.t_first_test.at(subj) << "\n";
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
    if(current_data_arg.dp_current.at(subj) < 10 || 
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

//[[Rcpp::export]]
double test_likelihood_calc(){
  viral_data_struct viral_data;
  current_data_struct current_data;
  current_parameters_struct current_parameters;
  settings_struct settings;
  
  viral_data.index = {0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5};
  viral_data.time = {-4,-3,-2,-1,0,-3,-2,-1,0,-2,-1,0,1,2,3,-3,-2,-1,0,1,2,-1,0,1,2,3,4,0,1,2,3,4};
  viral_data.viral_load = {0,8.2445,2.3393,18.2477,25.9738,0,12.9064,28.7029,27.1031,0,16.6276,29.0459,23.9215,25.7301,12.1522,12.6907,2.1197,18.6972,15.2154,9.7228,16.9616,39.806,22.8743,0,0,0,1.4154,20.4609,18.3101,12.7262,0,0};
  viral_data.subtype = {0,1,0,1,0,1};
  
  current_data.wp_current = {3,4,2,3,4,5};  
  current_data.dp_current = {35,34,33,20,40,30};
  current_data.tp_current = {1,0.5,0.1,-0.5,-1.5,-1};
  current_data.wr_current = {5,4,6,7,3,4};
  current_data.model_current = {1,1,2,2,3,3};
  
  current_parameters.wp_mean = {3.5,3.3};
  current_parameters.wp_sd = {1,1};
  current_parameters.tp_sd = {1.5,1.5};
  current_parameters.dp_mean = {36,34};
  current_parameters.dp_sd = {3,3};
  current_parameters.wr_mean = {5.5,5.3};
  current_parameters.wr_sd = {2,2};
  current_parameters.sigma = 5.0;
  
  settings.sensitivity = 0.99;
  
  double log_lh_test = log_likelihood(viral_data,
                                      current_data,
                                      current_parameters,
                                      settings);
  
  Rcout << "Total Likelihood: " << log_lh_test << "\n";
  
  return(log_lh_test);
}

//[[Rcpp::export]]
double test_likelihood_individual_calc(){
  viral_data_struct viral_data;
  current_data_struct current_data;
  current_parameters_struct current_parameters;
  settings_struct settings;
  
  viral_data.index = {0,0,0,0,0,
                      1,1,1,1,
                      2,2,2,2,2,2,
                      3,3,3,3,3,3,
                      4,4,4,4,4,4,
                      5,5,5,5,5};
  viral_data.time = {-4,-3,-2,-1,0,
                     -3,-2,-1,0,
                     -2,-1,0,1,2,3,
                     -3,-2,-1,0,1,2,
                     -1,0,1,2,3,4,
                     0,1,2,3,4};
  viral_data.viral_load = {0,8.2445,2.3393,18.2477,25.9738,
                           0,12.9064,28.7029,27.1031,
                           0,16.6276,29.0459,23.9215,25.7301,12.1522,
                           12.6907,2.1197,18.6972,15.2154,9.7228,16.9616,
                           39.806,22.8743,0,0,0,1.4154,
                           20.4609,18.3101,12.7262,0,0};
  viral_data.subtype = {0,1,0,1,0,1};
  
  current_data.wp_current = {3,4,2,3,4,5};  
  current_data.dp_current = {35,34,33,20,40,30};
  current_data.tp_current = {1,0.5,0.1,-0.5,-1.5,-1};
  current_data.wr_current = {5,4,6,7,3,4};
  current_data.model_current = {1,1,2,2,3,3};
  
  current_parameters.wp_mean = {3.5,3.3};
  current_parameters.wp_sd = {1,1};
  current_parameters.tp_sd = {1.5,1.5};
  current_parameters.dp_mean = {36,34};
  current_parameters.dp_sd = {3,3};
  current_parameters.wr_mean = {5.5,5.3};
  current_parameters.wr_sd = {2,2};
  current_parameters.sigma = 5.0;
  
  settings.sensitivity = 0.99;
  
  double log_lh_test = 0;
  
  for(int i = 0; i < 6; i++){
    log_lh_test += log_likelihood_subject(i,
                                          viral_data,
                                          current_data,
                                          current_parameters,
                                          settings);
    
    Rcout << "i " << i << " lh " << log_likelihood_subject(i,
                                                         viral_data,
                                                         current_data,
                                                         current_parameters,
                                                         settings) << "\n";
  }
  
  
  Rcout << "Total Lh: " << log_lh_test << "\n";
  
  return(log_lh_test);
}

//[[Rcpp::export]]
double log_likelihood_r_exp(std::vector<int> index_arg,
                            std::vector<double> viral_load_arg,
                            std::vector<double> time_arg,
                            std::vector<int> subtype_arg,
                            std::vector<double> t_first_positive_arg,
                            std::vector<double> t_last_positive_arg,
                            std::vector<double> t_first_test_arg,
                            std::vector<double> t_last_test_arg,
                            std::vector<int> n_positive_tests_arg,
                            std::vector<double> max_viral_load_arg,
                            double lod_arg,
                            double sensitivity_arg,
                            int n_subtypes_arg,
                            int n_subjects_arg,
                            int n_data_arg,
                            std::vector<double> wp_current_arg,
                            std::vector<double> tp_current_arg,
                            std::vector<double> dp_current_arg,
                            std::vector<double> wr_current_arg,
                            std::vector<int> model_current_arg,
                            double sigma_arg,
                            std::vector<double> wp_mean_arg,
                            std::vector<double> wp_sd_arg,
                            std::vector<double> tp_sd_arg,
                            std::vector<double> dp_mean_arg,
                            std::vector<double> dp_sd_arg,
                            std::vector<double> wr_mean_arg,
                            std::vector<double> wr_sd_arg
){
  // Organize data into structs
  viral_data_struct viral_data;
  current_data_struct current_data;
  current_parameters_struct current_parameters;
  settings_struct settings;
  
  make_viral_data_struct(viral_data,
                         index_arg,
                         viral_load_arg, 
                         time_arg, 
                         subtype_arg,
                         t_first_positive_arg,
                         t_last_positive_arg,
                         t_first_test_arg,
                         t_last_test_arg,
                         n_positive_tests_arg,
                         max_viral_load_arg);
  
  make_current_data_struct(current_data,
                           wp_current_arg,
                           tp_current_arg,
                           dp_current_arg,
                           wr_current_arg,
                           model_current_arg);
  
  make_current_parameters_struct(current_parameters,
                                 0.0,
                                 sigma_arg,
                                 wp_mean_arg,
                                 wp_sd_arg,
                                 tp_sd_arg,
                                 dp_mean_arg,
                                 dp_sd_arg,
                                 wr_mean_arg,
                                 wr_sd_arg);
  
  make_settings_struct(settings,
                       lod_arg,
                       sensitivity_arg,
                       0,
                       n_subtypes_arg,
                       n_subjects_arg,
                       n_data_arg);
  
  // Initialize likelihood
  current_parameters.log_likelihood = log_likelihood(viral_data,
                                                     current_data,
                                                     current_parameters,
                                                     settings);
  
  return(current_parameters.log_likelihood);
}