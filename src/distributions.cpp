#include <Rcpp.h>
#include "distributions.h"
using namespace Rcpp;

/* Random number generators */
double runif(double x0, double x1, rng_type& rng_arg) {
  boost::variate_generator<rng_type &, Dunif > rndm(rng_arg, Dunif(x0, x1));
  return rndm();
}

double rnorm_boost(double mean, double sd, rng_type& rng_arg) {
  boost::variate_generator<rng_type &, Dnorm > rndm(rng_arg, Dnorm(mean,sd));
  return rndm();
}

//[[Rcpp::export]]
void test_normal_pdf(){
  double value = 0.1;
  double mean = 1;
  double sd = 2;
  
  Rcout << "normal pdf mean=" << mean << ", sd=" << sd << ", x=" << value << ": " << pdf_normal(value,mean,sd) << "\n";
  
  value = 1.66;
  
  Rcout << "normal pdf mean=" << mean << ", sd=" << sd << ", x=" << value << ": " << pdf_normal(value,mean,sd) << "\n";
  
}

//[[Rcpp::export]]
void test_normal_cdf(){
  double value = 0.1;
  double mean = 1;
  double sd = 2;
  
  Rcout << "normal cdf mean=" << mean << ", sd=" << sd << ", x=" << value << ": " << cdf_normal(value,mean,sd) << "\n";
  
  value = 1.66;
  
  Rcout << "normal cdf mean=" << mean << ", sd=" << sd << ", x=" << value << ": " << cdf_normal(value,mean,sd) << "\n";
  
}

//[[Rcpp::export]]
void test_exp_pdf(){
  double value = 0.1;
  double rate = 0.1;
  
  Rcout << "exp pdf rate=" << rate << ", q=" << value << ": " << pdf_exponential(value,rate) << "\n";
  
  value = 0.001;
  
  Rcout << "exp pdf rate=" << rate << ", q=" << value << ": " << pdf_exponential(value,rate) << "\n";
  
}


//[[Rcpp::export]]
void test_exp_cdf(){
  double value = 0.1;
  double rate = 0.1;
  
  Rcout << "exp cdf rate=" << rate << ", q=" << value << ": " << cdf_exponential(value,rate) << "\n";
  
  value = 0.001;
  
  Rcout << "exp cdf rate=" << rate << ", q=" << value << ": " << cdf_exponential(value,rate) << "\n";
  
}

//[[Rcpp::export]]
void test_gamma_pdf(){
  double value = 4;
  double shape = 3;
  double scale = 2;
  
  Rcout << "gamma pdf shape=" << shape << " scale " << scale << ", q=" << value << ": " << pdf_gamma(value,shape,scale) << "\n";
  
  value = 2;
  
  Rcout << "gamma pdf shape=" << shape << " scale " << scale << ", q=" << value << ": " << pdf_gamma(value,shape,scale) << "\n";
  
}

//[[Rcpp::export]]
double test_rnorm_boost(double seed_arg){
  double mu = 0;
  double sigma = 1;
  
  rng_type rng_arg(seed_arg);
  
  return(rnorm_boost(mu, sigma, rng_arg));
}

//[[Rcpp::export]]
double test_runif(double seed_arg){
  double xmin = 1;
  double xmax = 3;
  
  rng_type rng_arg(seed_arg);
  
  return(runif(xmin,xmax,rng_arg));
}