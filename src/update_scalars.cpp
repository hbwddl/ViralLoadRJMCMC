#include <Rcpp.h>
#include "utilities.h"
#include "structs.h"
#include "distributions.h"
#include "log_likelihood.h"
#include "update_scalars.h"

using namespace Rcpp;

void update_wp_mean(int which_subtype_update,
                    current_data_struct& current_data_arg,
                    current_parameters_struct& current_parameters_arg,
                    settings_struct& settings_arg,
                    viral_data_struct& viral_data_arg,
                    priors_struct& priors_arg,
                    scaling_factors_struct& scaling_factors_arg,
                    double norm_0_1_draw,
                    double unif_0_1_draw){
  
  double wp_mean_proposed = current_parameters_arg.wp_mean.at(which_subtype_update) + (norm_0_1_draw*scaling_factors_arg.wp_mean_sf.at(which_subtype_update));
  
  if(wp_mean_proposed < 0 || wp_mean_proposed > priors_arg.wpmean_max){
    // out of bounds
    return;
  }
  
  current_parameters_struct parameters_proposed = current_parameters_arg;
  
  parameters_proposed.wp_mean.at(which_subtype_update) = wp_mean_proposed;
  
  double log_likelihood_proposed = current_parameters_arg.log_likelihood;
  
  log_likelihood_proposed = log_likelihood(viral_data_arg,
                                           current_data_arg,
                                           parameters_proposed,
                                           settings_arg);
  
  // Prior probability
  double log_prob_current = log(pdf_normal(current_parameters_arg.wp_mean.at(which_subtype_update),
                                       priors_arg.wpmean_mean,
                                       priors_arg.wpmean_sd));
  
  double log_prob_proposed = log(pdf_normal(wp_mean_proposed,
                                            priors_arg.wpmean_mean,
                                            priors_arg.wpmean_sd));
  
  // Acceptance probability
  double acp_pr = 0.0;
  
  acp_pr = exp(log_likelihood_proposed - current_parameters_arg.log_likelihood + log_prob_proposed - log_prob_current);
  
  if(unif_0_1_draw <= acp_pr){
    // Accept
    current_parameters_arg.log_likelihood = log_likelihood_proposed;
    current_parameters_arg.wp_mean.at(which_subtype_update) = wp_mean_proposed;
  } else{
    // Reject
    return;
  }
}