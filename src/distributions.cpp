#include <Rcpp.h>
#include "distributions.h"
using namespace Rcpp;

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