#include <Rcpp.h>
#include "utilities.h"
#include "structs.h"
#include "distributions.h"
#include "log_likelihood.h"
#include "update_scalars.h"

using namespace Rcpp;

#define gamma_prior_sd 0
#define uniform_prior_sd 1

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
    return;
  } else{
    // Reject
    return;
  }
}

void update_wp_sd(int which_subtype_update,
                    current_data_struct& current_data_arg,
                    current_parameters_struct& current_parameters_arg,
                    settings_struct& settings_arg,
                    viral_data_struct& viral_data_arg,
                    priors_struct& priors_arg,
                    scaling_factors_struct& scaling_factors_arg,
                    double norm_0_1_draw,
                    double unif_0_1_draw){
  double wp_sd_proposed = current_parameters_arg.wp_sd.at(which_subtype_update) + (norm_0_1_draw*scaling_factors_arg.wp_sd_sf.at(which_subtype_update));
  
  if(wp_sd_proposed < priors_arg.wpsd_min || wp_sd_proposed > priors_arg.wpsd_max){
    // out of bounds
    return;
  }
  
  current_parameters_struct parameters_proposed = current_parameters_arg;
  
  parameters_proposed.wp_sd.at(which_subtype_update) = wp_sd_proposed;
  
  double log_likelihood_proposed = current_parameters_arg.log_likelihood;
  
  log_likelihood_proposed = log_likelihood(viral_data_arg,
                                           current_data_arg,
                                           parameters_proposed,
                                           settings_arg);
  
  double log_prob_current = 0;
  double log_prob_proposed = 0;
  
  if(gamma_prior_sd == 1){
    log_prob_current = pdf_gamma(current_parameters_arg.wp_sd.at(which_subtype_update),priors_arg.wpsd_scale,1/priors_arg.wpsd_scale);
    log_prob_proposed = pdf_gamma(wp_sd_proposed,priors_arg.wpsd_scale,1/priors_arg.wpsd_scale);
  } else if(uniform_prior_sd == 1){
    log_prob_current = 0;
    log_prob_proposed = 0;
  } else{
    log_prob_current = pdf_exponential(current_parameters_arg.wp_sd.at(which_subtype_update),1/priors_arg.wpsd_scale);
    log_prob_proposed = pdf_exponential(wp_sd_proposed,1/priors_arg.wpsd_scale);
  }

  
  // Acceptance probability
  double acp_pr = 0.0;
  
  acp_pr = exp(log_likelihood_proposed - current_parameters_arg.log_likelihood + log_prob_proposed - log_prob_current);
  
  if(unif_0_1_draw <= acp_pr){
    // Accept
    current_parameters_arg.log_likelihood = log_likelihood_proposed;
    current_parameters_arg.wp_sd.at(which_subtype_update) = wp_sd_proposed;
    return;
  } else{
    // Reject
    return;
  }
}

void update_tp_sd(int which_subtype_update,
                  current_data_struct& current_data_arg,
                  current_parameters_struct& current_parameters_arg,
                  settings_struct& settings_arg,
                  viral_data_struct& viral_data_arg,
                  priors_struct& priors_arg,
                  scaling_factors_struct& scaling_factors_arg,
                  double norm_0_1_draw,
                  double unif_0_1_draw){
  double tp_sd_proposed = current_parameters_arg.tp_sd.at(which_subtype_update) + (norm_0_1_draw*scaling_factors_arg.tp_sd_sf.at(which_subtype_update));
  
  if(tp_sd_proposed < priors_arg.tpsd_min || tp_sd_proposed > priors_arg.tpsd_max){
    // out of bounds
    return;
  }
  
  current_parameters_struct parameters_proposed = current_parameters_arg;
  
  parameters_proposed.tp_sd.at(which_subtype_update) = tp_sd_proposed;
  
  double log_likelihood_proposed = current_parameters_arg.log_likelihood;
  
  log_likelihood_proposed = log_likelihood(viral_data_arg,
                                           current_data_arg,
                                           parameters_proposed,
                                           settings_arg);
  
  // Prior probability
  double log_prob_current = 0;
  double log_prob_proposed = 0;
  
  // Acceptance probability
  double acp_pr = 0.0;
  
  acp_pr = exp(log_likelihood_proposed - current_parameters_arg.log_likelihood + log_prob_proposed - log_prob_current);
  
  if(unif_0_1_draw <= acp_pr){
    // Accept
    current_parameters_arg.log_likelihood = log_likelihood_proposed;
    current_parameters_arg.tp_sd.at(which_subtype_update) = tp_sd_proposed;
    return;
  } else{
    // Reject
    return;
  }
}

void update_dp_mean(int which_subtype_update,
                  current_data_struct& current_data_arg,
                  current_parameters_struct& current_parameters_arg,
                  settings_struct& settings_arg,
                  viral_data_struct& viral_data_arg,
                  priors_struct& priors_arg,
                  scaling_factors_struct& scaling_factors_arg,
                  double norm_0_1_draw,
                  double unif_0_1_draw){
  double dp_mean_proposed = current_parameters_arg.dp_mean.at(which_subtype_update) + (norm_0_1_draw*scaling_factors_arg.dp_mean_sf.at(which_subtype_update));
  
  if(dp_mean_proposed < 0 || dp_mean_proposed > priors_arg.dpmean_max){
    // out of bounds
    return;
  }
  
  current_parameters_struct parameters_proposed = current_parameters_arg;
  
  parameters_proposed.dp_mean.at(which_subtype_update) = dp_mean_proposed;
  
  double log_likelihood_proposed = current_parameters_arg.log_likelihood;
  
  log_likelihood_proposed = log_likelihood(viral_data_arg,
                                           current_data_arg,
                                           parameters_proposed,
                                           settings_arg);
  
  // Prior probability
  double log_prob_current = log(pdf_normal(current_parameters_arg.dp_mean.at(which_subtype_update),
                                           priors_arg.dpmean_mean,
                                           priors_arg.dpmean_sd));
  
  double log_prob_proposed = log(pdf_normal(dp_mean_proposed,
                                            priors_arg.dpmean_mean,
                                            priors_arg.dpmean_sd));
  
  // Acceptance probability
  double acp_pr = 0.0;
  
  acp_pr = exp(log_likelihood_proposed - current_parameters_arg.log_likelihood + log_prob_proposed - log_prob_current);
  
  if(unif_0_1_draw <= acp_pr){
    // Accept
    current_parameters_arg.log_likelihood = log_likelihood_proposed;
    current_parameters_arg.dp_mean.at(which_subtype_update) = dp_mean_proposed;
    return;
  } else{
    // Reject
    return;
  }
}

void update_dp_sd(int which_subtype_update,
                  current_data_struct& current_data_arg,
                  current_parameters_struct& current_parameters_arg,
                  settings_struct& settings_arg,
                  viral_data_struct& viral_data_arg,
                  priors_struct& priors_arg,
                  scaling_factors_struct& scaling_factors_arg,
                  double norm_0_1_draw,
                  double unif_0_1_draw){
  double dp_sd_proposed = current_parameters_arg.dp_sd.at(which_subtype_update) + (norm_0_1_draw*scaling_factors_arg.dp_sd_sf.at(which_subtype_update));
  
  if(dp_sd_proposed < priors_arg.dpsd_min || dp_sd_proposed > priors_arg.dpsd_max){
    // out of bounds
    return;
  }
  
  current_parameters_struct parameters_proposed = current_parameters_arg;
  
  parameters_proposed.dp_sd.at(which_subtype_update) = dp_sd_proposed;
  
  double log_likelihood_proposed = current_parameters_arg.log_likelihood;
  
  log_likelihood_proposed = log_likelihood(viral_data_arg,
                                           current_data_arg,
                                           parameters_proposed,
                                           settings_arg);
  
  // Prior probability
  double log_prob_current = 0;
  double log_prob_proposed = 0;
  
  // Acceptance probability
  double acp_pr = 0.0;
  
  acp_pr = exp(log_likelihood_proposed - current_parameters_arg.log_likelihood + log_prob_proposed - log_prob_current);
  
  if(unif_0_1_draw <= acp_pr){
    // Accept
    current_parameters_arg.log_likelihood = log_likelihood_proposed;
    current_parameters_arg.dp_sd.at(which_subtype_update) = dp_sd_proposed;
    return;
  } else{
    // Reject
    return;
  }
}

void update_wr_mean(int which_subtype_update,
                  current_data_struct& current_data_arg,
                  current_parameters_struct& current_parameters_arg,
                  settings_struct& settings_arg,
                  viral_data_struct& viral_data_arg,
                  priors_struct& priors_arg,
                  scaling_factors_struct& scaling_factors_arg,
                  double norm_0_1_draw,
                  double unif_0_1_draw){
  double wr_mean_proposed = current_parameters_arg.wr_mean.at(which_subtype_update) + (norm_0_1_draw*scaling_factors_arg.wr_mean_sf.at(which_subtype_update));
  
  if(wr_mean_proposed < 0 || wr_mean_proposed > priors_arg.wrmean_max){
    // out of bounds
    return;
  }
  
  current_parameters_struct parameters_proposed = current_parameters_arg;
  
  parameters_proposed.wr_mean.at(which_subtype_update) = wr_mean_proposed;
  
  double log_likelihood_proposed = current_parameters_arg.log_likelihood;
  
  log_likelihood_proposed = log_likelihood(viral_data_arg,
                                           current_data_arg,
                                           parameters_proposed,
                                           settings_arg);
  
  // Prior probability
  double log_prob_current = log(pdf_normal(current_parameters_arg.wr_mean.at(which_subtype_update),
                                           priors_arg.wrmean_mean,
                                           priors_arg.wrmean_sd));
  
  double log_prob_proposed = log(pdf_normal(wr_mean_proposed,
                                            priors_arg.wrmean_mean,
                                            priors_arg.wrmean_sd));
  
  // Acceptance probability
  double acp_pr = 0.0;
  
  acp_pr = exp(log_likelihood_proposed - current_parameters_arg.log_likelihood + log_prob_proposed - log_prob_current);
  
  if(unif_0_1_draw <= acp_pr){
    // Accept
    current_parameters_arg.log_likelihood = log_likelihood_proposed;
    current_parameters_arg.wr_mean.at(which_subtype_update) = wr_mean_proposed;
    return;
  } else{
    // Reject
    return;
  }
}

void update_wr_sd(int which_subtype_update,
                    current_data_struct& current_data_arg,
                    current_parameters_struct& current_parameters_arg,
                    settings_struct& settings_arg,
                    viral_data_struct& viral_data_arg,
                    priors_struct& priors_arg,
                    scaling_factors_struct& scaling_factors_arg,
                    double norm_0_1_draw,
                    double unif_0_1_draw){
  double wr_sd_proposed = current_parameters_arg.wr_sd.at(which_subtype_update) + (norm_0_1_draw*scaling_factors_arg.wr_sd_sf.at(which_subtype_update));
  
  if(wr_sd_proposed < priors_arg.wrsd_min || wr_sd_proposed > priors_arg.wrsd_max){
    // out of bounds
    return;
  }
  
  current_parameters_struct parameters_proposed = current_parameters_arg;
  
  parameters_proposed.wr_sd.at(which_subtype_update) = wr_sd_proposed;
  
  double log_likelihood_proposed = current_parameters_arg.log_likelihood;
  
  log_likelihood_proposed = log_likelihood(viral_data_arg,
                                           current_data_arg,
                                           parameters_proposed,
                                           settings_arg);
  
  // Prior probability
  
  double log_prob_current = 0;
  double log_prob_proposed = 0;
  
  if(gamma_prior_sd == 1){
    log_prob_current = pdf_gamma(current_parameters_arg.wr_sd.at(which_subtype_update),priors_arg.wrsd_scale,1/priors_arg.wrsd_scale);
    log_prob_proposed = pdf_gamma(wr_sd_proposed,priors_arg.wrsd_scale,1/priors_arg.wrsd_scale);
  } else if(uniform_prior_sd == 1){
    log_prob_current = 0;
    log_prob_proposed = 0;
  } else{
    log_prob_current = pdf_exponential(current_parameters_arg.wr_sd.at(which_subtype_update),1/priors_arg.wrsd_scale);
    log_prob_proposed = pdf_exponential(wr_sd_proposed,1/priors_arg.wrsd_scale);
  }
  
  // Acceptance probability
  double acp_pr = 0.0;
  
  acp_pr = exp(log_likelihood_proposed - current_parameters_arg.log_likelihood + log_prob_proposed - log_prob_current);
  
  if(unif_0_1_draw <= acp_pr){
    // Accept
    current_parameters_arg.log_likelihood = log_likelihood_proposed;
    current_parameters_arg.wr_sd.at(which_subtype_update) = wr_sd_proposed;
    return;
  } else{
    // Reject
    return;
  }
}

void update_sigma(current_data_struct& current_data_arg,
                  current_parameters_struct& current_parameters_arg,
                  settings_struct& settings_arg,
                  viral_data_struct& viral_data_arg,
                  priors_struct& priors_arg,
                  scaling_factors_struct& scaling_factors_arg,
                  double norm_0_1_draw,
                  double unif_0_1_draw){
  double sigma_proposed = current_parameters_arg.sigma + (norm_0_1_draw*scaling_factors_arg.sigma_sf);
  
  if(sigma_proposed < priors_arg.sigma_min || sigma_proposed > priors_arg.sigma_max){
    // out of bounds
    return;
  }
  
  current_parameters_struct parameters_proposed = current_parameters_arg;
  
  parameters_proposed.sigma = sigma_proposed;
  
  double log_likelihood_proposed = current_parameters_arg.log_likelihood;
  
  log_likelihood_proposed = log_likelihood(viral_data_arg,
                                           current_data_arg,
                                           parameters_proposed,
                                           settings_arg);
  
  // Prior probability
  double log_prob_current = 0;
  double log_prob_proposed = 0;
  
  // Acceptance probability
  double acp_pr = 0.0;
  
  acp_pr = exp(log_likelihood_proposed - current_parameters_arg.log_likelihood + log_prob_proposed - log_prob_current);
  
  if(unif_0_1_draw <= acp_pr){
    // Accept
    current_parameters_arg.log_likelihood = log_likelihood_proposed;
    current_parameters_arg.sigma = sigma_proposed;
    return;
  } else{
    // Reject
    return;
  }
}