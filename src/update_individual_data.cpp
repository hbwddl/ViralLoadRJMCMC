#include <Rcpp.h>
#include "structs.h"
#include "utilities.h"
#include "distributions.h"
#include "log_likelihood.h"
#include "debug.h"
#include "update_individual_data.h"

using namespace Rcpp;

void update_wp_i(int index_update,
                 current_data_struct& current_data_arg,
                 current_parameters_struct& current_parameters_arg,
                 settings_struct& settings_arg,
                 viral_data_struct& viral_data_arg,
                 priors_struct& priors_arg,
                 scaling_factors_struct& scaling_factors_arg,
                 double norm_0_1_draw,
                 double unif_0_1_draw){
  int model_current_i = current_data_arg.model_current.at(index_update);
  if(model_current_i == 3){
    return;
  }
  
  double wp_current_i = current_data_arg.wp_current.at(index_update);
  double tp_current_i = current_data_arg.tp_current.at(index_update);
  double first_gt0_i = viral_data_arg.t_first_positive.at(index_update);
  double last_gt0_i = viral_data_arg.t_last_positive.at(index_update);
  int subtype_i = viral_data_arg.subtype.at(index_update);
  
  double wp_proposed_i = wp_current_i + norm_0_1_draw*scaling_factors_arg.wp_sf;
  
  double log_lh_proposed = current_parameters_arg.log_likelihood;
  
  current_data_struct proposed_data = current_data_arg;
  
  proposed_data.wp_current.at(index_update) = wp_proposed_i;
  
  double wp_min_i = priors_arg.wp_min;
  double wp_max_i = priors_arg.wp_max;
  
  if(wp_current_i < priors_arg.wp_min){
    wp_proposed_i = priors_arg.wp_min + 0.1;
    proposed_data.wp_current.at(index_update) = wp_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.wp_current.at(index_update) = wp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  }
  
  if(wp_current_i > priors_arg.wp_max){
    wp_proposed_i = priors_arg.wp_max - 0.1;
    proposed_data.wp_current.at(index_update) = wp_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.wp_current.at(index_update) = wp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  }
  
  if(wp_proposed_i > priors_arg.wp_max || wp_proposed_i < priors_arg.wp_min){
    return;
  }
  
  log_lh_proposed = log_likelihood(viral_data_arg,
                                   proposed_data,
                                   current_parameters_arg,
                                   settings_arg);
  
  double acp_pr = exp(log_lh_proposed - current_parameters_arg.log_likelihood);
  
  if(acp_pr < unif_0_1_draw){
    // Accept
    current_data_arg.wp_current.at(index_update) = wp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    
    return;
  } else{
    return;
  }
  
}

void update_tp_i(int index_update,
                 current_data_struct& current_data_arg,
                 current_parameters_struct& current_parameters_arg,
                 settings_struct& settings_arg,
                 viral_data_struct& viral_data_arg,
                 priors_struct& priors_arg,
                 scaling_factors_struct& scaling_factors_arg,
                 double norm_0_1_draw,
                 double unif_0_1_draw){
  int model_current_i = current_data_arg.model_current.at(index_update);
  
  double tp_min_i = -priors_arg.wp_max;
  double tp_max_i = priors_arg.wr_max;
  
  double wp_current_i = current_data_arg.wp_current.at(index_update);
  double tp_current_i = current_data_arg.tp_current.at(index_update);
  double first_gt0_i = viral_data_arg.t_first_positive.at(index_update);
  double last_gt0_i = viral_data_arg.t_last_positive.at(index_update);
  int subtype_i = viral_data_arg.subtype.at(index_update);
  
  double tp_proposed_i = tp_current_i + norm_0_1_draw*scaling_factors_arg.tp_sf;
  
  double log_lh_proposed = current_parameters_arg.log_likelihood;
  
  current_data_struct proposed_data = current_data_arg;
  
  proposed_data.tp_current.at(index_update) = tp_proposed_i;
  
  if(model_current_i == 1){ // Increase only, tp_min is last test
    tp_min_i = viral_data_arg.t_last_test.at(index_update);
  } else if(model_current_i == 2){ // Peak, tp bounded by tests
    tp_min_i = viral_data_arg.t_first_test.at(index_update);
    tp_max_i = viral_data_arg.t_last_test.at(index_update);
  } else if(model_current_i == 3){ // Only decline, tp bounded above by first test
    tp_max_i = viral_data_arg.t_first_test.at(index_update);
  } else{
    Rcout << "ERR ";
    print_pos(__FILE__,
              __LINE__,
              1);
  }
  
  if(tp_current_i < tp_min_i){
    tp_proposed_i = tp_min_i + 0.1;
    proposed_data.tp_current.at(index_update) = tp_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.tp_current.at(index_update) = tp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  }
  
  if(tp_current_i > tp_max_i){
    tp_proposed_i = tp_max_i - 0.1;
    proposed_data.tp_current.at(index_update) = tp_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.tp_current.at(index_update) = tp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  }
  
  if(tp_proposed_i > tp_max_i || tp_proposed_i < tp_min_i){
    return;
  }
  
  log_lh_proposed = log_likelihood(viral_data_arg,
                                   proposed_data,
                                   current_parameters_arg,
                                   settings_arg);
  
  double acp_pr = exp(log_lh_proposed - current_parameters_arg.log_likelihood);
  
  if(acp_pr < unif_0_1_draw){
    // Accept
    current_data_arg.tp_current.at(index_update) = tp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    
    return;
  } else{
    return;
  }
  
}

void update_dp_i(int index_update,
                 current_data_struct& current_data_arg,
                 current_parameters_struct& current_parameters_arg,
                 settings_struct& settings_arg,
                 viral_data_struct& viral_data_arg,
                 priors_struct& priors_arg,
                 scaling_factors_struct& scaling_factors_arg,
                 double norm_0_1_draw,
                 double unif_0_1_draw){
  
  double dp_min_i = viral_data_arg.max_viral_load.at(index_update)/2;
  
  double dp_current_i = current_data_arg.dp_current.at(index_update);
  
  double dp_proposed_i = current_data_arg.dp_current.at(index_update) + norm_0_1_draw*scaling_factors_arg.dp_sf;
  
  double log_lh_proposed = current_parameters_arg.log_likelihood;
  
  current_data_struct proposed_data = current_data_arg;
  
  if(dp_current_i > settings_arg.lod){
    dp_proposed_i = settings_arg.lod - 0.1;
    proposed_data.dp_current.at(index_update) = dp_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.dp_current.at(index_update) = dp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  }
  
  if(dp_current_i < dp_min_i){
    dp_proposed_i = dp_min_i + 0.1;
    proposed_data.dp_current.at(index_update) = dp_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.dp_current.at(index_update) = dp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  }
  
  
  if(dp_proposed_i > settings_arg.lod || dp_proposed_i < dp_min_i){
    return;
  }
  
  proposed_data.dp_current.at(index_update) = dp_proposed_i;

  log_lh_proposed = log_likelihood(viral_data_arg,
                                   proposed_data,
                                   current_parameters_arg,
                                   settings_arg);
  
  double acp_pr = exp(log_lh_proposed - current_parameters_arg.log_likelihood);
  
  if(acp_pr < unif_0_1_draw){
    // Accept
    current_data_arg.dp_current.at(index_update) = dp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    
    return;
  } else{
    return;
  }
  
}

void update_wr_i(int index_update,
                 current_data_struct& current_data_arg,
                 current_parameters_struct& current_parameters_arg,
                 settings_struct& settings_arg,
                 viral_data_struct& viral_data_arg,
                 priors_struct& priors_arg,
                 scaling_factors_struct& scaling_factors_arg,
                 double norm_0_1_draw,
                 double unif_0_1_draw){
  
  int model_current_i = current_data_arg.model_current.at(index_update);
  
  if(model_current_i == 1){ // Don't update if only proliferation stage
    return;
  }
  
  double wr_current_i = current_data_arg.wr_current.at(index_update);
  double tp_current_i = current_data_arg.tp_current.at(index_update);
  double first_gt0_i = viral_data_arg.t_first_positive.at(index_update);
  double last_gt0_i = viral_data_arg.t_last_positive.at(index_update);
  int subtype_i = viral_data_arg.subtype.at(index_update);
  
  double wr_proposed_i = wr_current_i + norm_0_1_draw*scaling_factors_arg.wr_sf;
  
  double log_lh_proposed = current_parameters_arg.log_likelihood;
  
  current_data_struct proposed_data = current_data_arg;
  
  proposed_data.wr_current.at(index_update) = wr_proposed_i;
  
  if(wr_current_i < priors_arg.wr_min){
    wr_proposed_i = priors_arg.wr_min + 0.1;
    proposed_data.wr_current.at(index_update) = wr_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.wr_current.at(index_update) = wr_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  }
  
  if(wr_current_i > priors_arg.wr_max){
    wr_proposed_i = priors_arg.wr_max - 0.1;
    proposed_data.wr_current.at(index_update) = wr_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.wr_current.at(index_update) = wr_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  }
  
  if(wr_proposed_i > priors_arg.wr_max || wr_proposed_i < priors_arg.wr_min){
    return;
  }
  
  log_lh_proposed = log_likelihood(viral_data_arg,
                                   proposed_data,
                                   current_parameters_arg,
                                   settings_arg);
  
  double acp_pr = exp(log_lh_proposed - current_parameters_arg.log_likelihood);
  
  if(acp_pr < unif_0_1_draw){
    // Accept
    current_data_arg.wr_current.at(index_update) = wr_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    
    return;
  } else{
    return;
  }
  
}