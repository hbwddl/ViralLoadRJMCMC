#include <Rcpp.h>
#include "structs.h"
#include "utilities.h"
#include "distributions.h"
#include "log_likelihood.h"
#include "debug.h"
#include "update_individual_data.h"

using namespace Rcpp;

#define debug_wp_update 0
#define debug_tp_update 0
#define debug_dp_update 0
#define debug_wr_update 0
#define debug_model_update 0

#define tp_unbounded 0

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
  
  if(debug_wp_update == 1){
    Rcout << "index " << index_update << " wp_current_i " << wp_current_i << " norm_0_1_draw " << norm_0_1_draw << " sf " << scaling_factors_arg.wp_sf << " wp_proposed_i " << wp_proposed_i << " ";
  }

  double log_lh_proposed = current_parameters_arg.log_likelihood;
  
  current_data_struct proposed_data = current_data_arg;
  
  proposed_data.wp_current.at(index_update) = wp_proposed_i;
  
  double wp_min_i = priors_arg.wp_min;
  double wp_max_i = priors_arg.wp_max;

  wp_min_i = std::max(priors_arg.wp_min,tp_current_i-first_gt0_i);
  wp_max_i = std::min(priors_arg.wp_max,(tp_current_i-first_gt0_i)+2);

  if(debug_wp_update == 1){
    Rcout << "wp_min_i " << wp_min_i << " wp_max_i " << wp_max_i << " ";
  }
  
  if(wp_max_i < wp_min_i){
    wp_min_i = priors_arg.wp_min;
    wp_max_i = priors_arg.wp_max;
  }
  
  if(wp_current_i < wp_min_i){
    wp_proposed_i = wp_min_i + 0.0001;
    proposed_data.wp_current.at(index_update) = wp_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.wp_current.at(index_update) = wp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    
    if(debug_wp_update == 1){
      Rcout << " wp current lo fix " << wp_proposed_i << " " << current_data_arg.wp_current.at(index_update) << "\n";
    }
    
    return;
  }
  
  if(wp_current_i > wp_max_i){
    wp_proposed_i = wp_max_i - 0.0001;
    proposed_data.wp_current.at(index_update) = wp_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.wp_current.at(index_update) = wp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    
    if(debug_wp_update == 1){
      Rcout << " wp current hi fix\n";
    }
    
    return;
  }

  // if((tp_current_i - wp_current_i) > first_gt0_i){
  //   wp_proposed_i = (tp_current_i - first_gt0_i) + 0.0001;
  // 
  //   if(wp_proposed_i < 0){
  //     Rcout << "ERR " << __FILE__ << " " << __LINE__ << "\n";
  //   }
  // 
  //   proposed_data.wp_current.at(index_update) = wp_proposed_i;
  // 
  //   log_lh_proposed = log_likelihood(viral_data_arg,
  //                                    proposed_data,
  //                                    current_parameters_arg,
  //                                    settings_arg);
  // 
  //   current_data_arg.wp_current.at(index_update) = wp_proposed_i;
  //   current_parameters_arg.log_likelihood = log_lh_proposed;
  // 
  //   if(debug_wp_update == 1){
  //     Rcout << " wp current span fix\n";
  //   }
  // 
  //   return;
  // }
  
  if(wp_proposed_i > wp_max_i || wp_proposed_i < wp_min_i){
    if(debug_wp_update == 1){
      Rcout << " invalid proposal\n";
    }
    return;
  }
  
  log_lh_proposed = current_parameters_arg.log_likelihood -
                    log_likelihood_subject(index_update,
                                           viral_data_arg,
                                           current_data_arg,
                                           current_parameters_arg,
                                           settings_arg) +
                     log_likelihood_subject(index_update,
                                            viral_data_arg,
                                            proposed_data,
                                            current_parameters_arg,
                                            settings_arg);
  
  if(debug_wp_update == 1){
    Rcout << " current_parameters_arg.wp_mean " << current_parameters_arg.wp_mean.at(subtype_i) << " current_parameters_arg.wp_sd " << current_parameters_arg.wp_sd.at(subtype_i) << " log_lh_proposed " << log_lh_proposed << " current_parameters_arg.log_likelihood " << current_parameters_arg.log_likelihood << " ";
  }
  
  double acp_pr = exp(log_lh_proposed - current_parameters_arg.log_likelihood);
  
  if(debug_wp_update == 1){
    Rcout << " acp_pr " << acp_pr << " unif_0_1_draw " << unif_0_1_draw << " ";
  }
  
  if(unif_0_1_draw < acp_pr){
    // Accept
    current_data_arg.wp_current.at(index_update) = wp_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    
    if(debug_wp_update == 1){
      Rcout << "\n";
    }
    
    return;
  } else{
    if(debug_wp_update == 1){
      Rcout << "\n";
    }
    
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
  
  // double tp_min_i = -2.0;
  // double tp_max_i = 2.0;
  
  double wp_current_i = current_data_arg.wp_current.at(index_update);
  double wr_current_i = current_data_arg.wr_current.at(index_update);
  double tp_current_i = current_data_arg.tp_current.at(index_update);
  double first_gt0_i = viral_data_arg.t_first_positive.at(index_update);
  double last_gt0_i = viral_data_arg.t_last_positive.at(index_update);
  double first_test_i = viral_data_arg.t_first_test.at(index_update);
  double last_test_i = viral_data_arg.t_last_test.at(index_update);
  int subtype_i = viral_data_arg.subtype.at(index_update);
  
  double tp_proposed_i = tp_current_i + norm_0_1_draw*scaling_factors_arg.tp_sf;
  // normal draw
  tp_proposed_i = norm_0_1_draw*current_parameters_arg.tp_sd.at(0);
  
  double log_lh_proposed = current_parameters_arg.log_likelihood;
  
  current_data_struct proposed_data = current_data_arg;
  
  proposed_data.tp_current.at(index_update) = tp_proposed_i;
  
  if(model_current_i == 1){ // Proliferation only, tp_min is last test
    tp_min_i = last_test_i;
    tp_max_i = last_test_i + 2;
    // tp_max_i = first_test_i + priors_arg.wp_max;
  } else if(model_current_i == 2){ // Peak, tp bounded by tests
    // tp_min_i = std::max(std::max(first_gt0_i,-2.0),last_gt0_i - wr_current_i);
    // tp_max_i = std::min(std::min(last_gt0_i,2.0),first_gt0_i + wp_current_i);

    tp_min_i = first_test_i;
    tp_max_i = last_test_i;
    
  } else if(model_current_i == 3){ // Only clearance, tp bounded above by first test
    tp_max_i = first_test_i;
    tp_min_i = first_test_i - 2;
    // tp_min_i = first_test_i - priors_arg.wr_max;
  } else{
    Rcout << "ERR ";
    print_pos(__FILE__,
              __LINE__,
              1);
  }
  
  if(tp_unbounded == 1){
    tp_min_i = -priors_arg.wp_max;
    tp_max_i = priors_arg.wr_max;
  }
  
  if(tp_current_i < tp_min_i){
    tp_proposed_i = tp_min_i + 0.0001;
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
    tp_proposed_i = tp_max_i - 0.0001;
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
  
  log_lh_proposed = current_parameters_arg.log_likelihood - 
                              log_likelihood_subject(index_update,
                                                     viral_data_arg,
                                                     current_data_arg,
                                                     current_parameters_arg,
                                                     settings_arg) + 
                             log_likelihood_subject(index_update,
                                                    viral_data_arg,
                                                    proposed_data,
                                                    current_parameters_arg,
                                                    settings_arg);
  
  double acp_pr = exp(log_lh_proposed - current_parameters_arg.log_likelihood);
  
  if(unif_0_1_draw < acp_pr){
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
  double dp_min_i = 0.0;
  
  // dp_min_i = settings_arg.lod*(1.0/2.0);
  
  // dp_min_i = viral_data_arg.max_viral_load.at(index_update) * (1.0/2.0);
  
  // dp_min_i = viral_data_arg.max_viral_load.at(index_update) - 3*current_parameters_arg.dp_sd.at(0);
  
  double dp_current_i = current_data_arg.dp_current.at(index_update);
  
  double dp_proposed_i = current_data_arg.dp_current.at(index_update) + norm_0_1_draw*scaling_factors_arg.dp_sf;
  
  double log_lh_proposed = current_parameters_arg.log_likelihood;
  
  current_data_struct proposed_data = current_data_arg;
  
  if(dp_current_i > settings_arg.lod){
    dp_proposed_i = settings_arg.lod - 0.0001;
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
    dp_proposed_i = dp_min_i + 0.0001;
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

  log_lh_proposed = current_parameters_arg.log_likelihood - 
    log_likelihood_subject(index_update,
                           viral_data_arg,
                           current_data_arg,
                           current_parameters_arg,
                           settings_arg) + 
                             log_likelihood_subject(index_update,
                                                    viral_data_arg,
                                                    proposed_data,
                                                    current_parameters_arg,
                                                    settings_arg);
  
  double acp_pr = exp(log_lh_proposed - current_parameters_arg.log_likelihood);
  
  if(unif_0_1_draw < acp_pr){
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
  
  double wr_min_i = priors_arg.wr_min;
  double wr_max_i = priors_arg.wr_max;
  
  wr_min_i = std::max(priors_arg.wr_min,last_gt0_i-tp_current_i);
  // wr_max_i = std::min(priors_arg.wr_max,(last_gt0_i-tp_current_i)+2);
  // // 
  
  if(debug_wr_update == 1){
    Rcout << "index " << index_update << " wr_current_i " << wr_current_i << " norm_0_1_draw " << norm_0_1_draw << " sf " << scaling_factors_arg.wr_sf << " wr_proposed_i " << wr_proposed_i << " ";
  }
  
  if(debug_wp_update == 1){
    Rcout << "wr_min_i " << wr_min_i << " wr_max_i " << wr_max_i << " ";
  }
  
  if(wr_max_i < wr_min_i){
    wr_min_i = priors_arg.wr_min;
    wr_max_i = priors_arg.wr_max;
  }
  
  if(wr_current_i < wr_min_i){
    wr_proposed_i = wr_min_i + 0.0001;
    proposed_data.wr_current.at(index_update) = wr_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.wr_current.at(index_update) = wr_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  }
  
  if(wr_current_i > wr_max_i){
    wr_proposed_i = wr_max_i - 0.0001;
    proposed_data.wr_current.at(index_update) = wr_proposed_i;
    
    log_lh_proposed = log_likelihood(viral_data_arg,
                                     proposed_data,
                                     current_parameters_arg,
                                     settings_arg);
    
    current_data_arg.wr_current.at(index_update) = wr_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  }
  
  // if(tp_current_i + wr_current_i < last_gt0_i){
  //   wr_proposed_i = (last_gt0_i - tp_current_i) + 0.01;
  // 
  //   proposed_data.wr_current.at(index_update) = wr_proposed_i;
  // 
  //   log_lh_proposed = log_likelihood(viral_data_arg,
  //                                    proposed_data,
  //                                    current_parameters_arg,
  //                                    settings_arg);
  // 
  //   current_data_arg.wr_current.at(index_update) = wr_proposed_i;
  //   current_parameters_arg.log_likelihood = log_lh_proposed;
  //   return;
  // }

  
  if(wr_proposed_i > wr_max_i || wr_proposed_i < wr_min_i){
    return;
  }
  
  log_lh_proposed = current_parameters_arg.log_likelihood - 
    log_likelihood_subject(index_update,
                           viral_data_arg,
                           current_data_arg,
                           current_parameters_arg,
                           settings_arg) + 
                             log_likelihood_subject(index_update,
                                                    viral_data_arg,
                                                    proposed_data,
                                                    current_parameters_arg,
                                                    settings_arg);
  
  double acp_pr = exp(log_lh_proposed - current_parameters_arg.log_likelihood);
  
  
  if(debug_wr_update == 1){
    Rcout << "\n";
  }
  
  if(unif_0_1_draw < acp_pr){
    // Accept
    current_data_arg.wr_current.at(index_update) = wr_proposed_i;
    current_parameters_arg.log_likelihood = log_lh_proposed;
    
    return;
  } else{
    return;
  }
  
}