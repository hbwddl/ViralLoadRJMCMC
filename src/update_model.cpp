#include <Rcpp.h>
#include "structs.h"
#include "utilities.h"
#include "distributions.h"
#include "log_likelihood.h"
#include "debug.h"
#include "update_individual_data.h"

using namespace Rcpp;

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
                 double unif_0_1_draw_acp){
  int model_current_i = current_data_arg.model_current.at(index_update);
  double wp_current_i = current_data_arg.wp_current.at(index_update);
  double tp_current_i = current_data_arg.tp_current.at(index_update);
  double dp_current_i = current_data_arg.dp_current.at(index_update);
  double wr_current_i = current_data_arg.wr_current.at(index_update);
  double first_test_i = viral_data_arg.t_first_test.at(index_update);
  double last_test_i = viral_data_arg.t_last_test.at(index_update);
  double first_gt0_i = viral_data_arg.t_first_positive.at(index_update);
  double last_gt0_i = viral_data_arg.t_last_positive.at(index_update);
  int subtype_i = viral_data_arg.subtype.at(index_update);
  
  double log_lh_proposed = current_parameters_arg.log_likelihood;
  double acp_pr_multiply = 1;
  double acp_pr = 0;
  current_data_struct proposed_data = current_data_arg;
  
  int model_proposed_i = model_current_i;
  
  if(unif_0_1_draw_model < 0.3333 && unif_0_1_draw_model >= 0){
    /* Selecting model 1 */
    model_proposed_i = 1;
  } else if(unif_0_1_draw_model < 0.6667 && unif_0_1_draw_model >= 0.3333){
    /* Selecting model 2 */
    model_proposed_i = 2;
  } else if(unif_0_1_draw_model <= 1 && unif_0_1_draw_model >= 0.6667){
    /* Selecting model 3 */
    model_proposed_i = 3;
  } else{
    print_pos(__FILE__,
              __LINE__,
              1);
  }
  
  if(model_current_i == model_proposed_i){
    return;
  }
  
  double wp_proposed_i = wp_current_i;
  double tp_proposed_i = tp_current_i;
  double dp_proposed_i = dp_current_i;
  double wr_proposed_i = wr_current_i;
  double u_proposed = 0.5;
  double v_proposed = 0.5;
  
  if(model_current_i == 1 && model_proposed_i == 2){
    wp_proposed_i = wp_current_i - (2*(tp_current_i - last_test_i));
    tp_proposed_i = 2*last_test_i - tp_current_i;
    dp_proposed_i = dp_current_i;
    
    u_proposed = unif_0_1_draw_rjmcmc_u;
    wr_proposed_i = 2*(tp_current_i - last_test_i) - log(u_proposed);
    
    acp_pr_multiply = 1/u_proposed;
    
  } else if(model_current_i == 2 && model_proposed_i == 1){
    wp_proposed_i = wp_current_i + (2*(last_test_i - tp_current_i));
    tp_proposed_i = 2*last_test_i - tp_current_i;
    dp_proposed_i = dp_current_i;
    u_proposed = exp(2*(last_test_i - tp_current_i) - wr_current_i);
    
    acp_pr_multiply = exp(2*(last_test_i - tp_current_i) - wr_current_i);
    
  } else if(model_current_i == 2 && model_proposed_i == 3){
    v_proposed = exp(2*(tp_current_i - first_test_i) - wp_current_i);
    tp_proposed_i = 2*first_test_i - tp_current_i;
    dp_proposed_i = dp_current_i;
    wr_proposed_i = wr_current_i + 2*(tp_current_i - first_test_i);
    
    acp_pr_multiply = exp(2*(tp_current_i - first_test_i) - wp_current_i);
    
  } else if(model_current_i == 3 && model_proposed_i == 2){
    v_proposed = unif_0_1_draw_rjmcmc_v;
    wp_proposed_i = 2*(first_test_i - tp_current_i) - log(v_proposed);
    tp_proposed_i = 2*first_test_i - tp_current_i;
    dp_proposed_i = dp_current_i;
    wr_proposed_i = wr_current_i + 2*(tp_current_i - first_test_i);
    
    acp_pr_multiply = 1/v_proposed;
    
  } else if(model_current_i == 1 && model_proposed_i == 3){
    u_proposed = unif_0_1_draw_rjmcmc_u;
    v_proposed = exp((last_test_i - first_test_i) - wp_current_i);
    tp_proposed_i = first_test_i + last_test_i - tp_current_i;
    dp_proposed_i = dp_current_i;
    wr_proposed_i = (last_test_i - first_test_i) - log(u_proposed);
    
    acp_pr_multiply = exp((last_test_i - first_test_i) - wp_current_i)/u_proposed;
    
  } else if(model_current_i == 3 && model_proposed_i == 1){
    u_proposed = exp((last_test_i - first_test_i) - wr_current_i);
    v_proposed = unif_0_1_draw_rjmcmc_v;
    wp_proposed_i = (last_test_i - first_test_i) - log(v_proposed);
    tp_proposed_i = first_test_i + last_test_i - tp_current_i;
    dp_proposed_i = dp_current_i;
    
    acp_pr_multiply = exp((last_test_i - first_test_i) - wr_current_i) / v_proposed;
    
  } else{
    print_pos(__FILE__,
              __LINE__,
              1);
    
    Rcout << "Model current " << model_current_i << " model proposed " << model_proposed_i << "\n";
  }
  
  proposed_data.model_current.at(index_update) = model_proposed_i;
  proposed_data.wp_current.at(index_update) = wp_proposed_i;
  proposed_data.tp_current.at(index_update) = tp_proposed_i;
  proposed_data.dp_current.at(index_update) = dp_proposed_i;
  proposed_data.wr_current.at(index_update) = wr_proposed_i;
  
  int data_err = check_data_err(viral_data_arg,
                                proposed_data,
                                current_parameters_arg,
                                settings_arg,
                                priors_arg);
  
  if(data_err == 1){
    return;
  }
  
  log_lh_proposed = log_likelihood(viral_data_arg,
                                   proposed_data,
                                   current_parameters_arg,
                                   settings_arg);
  
  acp_pr = exp(log_lh_proposed - current_parameters_arg.log_likelihood)*acp_pr_multiply;
  
  if(unif_0_1_draw_acp < acp_pr){
    // Accept
    current_data_arg.model_current.at(index_update) = proposed_data.model_current.at(index_update);
    current_data_arg.wp_current.at(index_update) = proposed_data.wp_current.at(index_update);
    current_data_arg.tp_current.at(index_update) = proposed_data.tp_current.at(index_update);
    current_data_arg.dp_current.at(index_update) = proposed_data.dp_current.at(index_update);
    current_data_arg.wr_current.at(index_update) = proposed_data.wr_current.at(index_update);
    
    current_parameters_arg.log_likelihood = log_lh_proposed;
    return;
  } else{
    return;
  }
  
}
