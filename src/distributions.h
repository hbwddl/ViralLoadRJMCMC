#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <boost/math/distributions.hpp>

typedef boost::math::normal_distribution<double> normal_mdist; //gamma_mdist(shape, scale)
typedef boost::math::exponential_distribution<double> exp_mdist; //exp_mdist(rate)

/* Normal PDF */
inline double pdf_normal(double x, double mean, double sd){
  double pdf_normal = 0.0;
  
  pdf_normal = pdf(normal_mdist(mean,sd),x);
  
  if (pdf_normal == 0){
    double dbl_min = std::numeric_limits< double >::min();
    pdf_normal = dbl_min;
  }
  
  return(pdf_normal);
}

inline double cdf_normal(double q, double mean, double sd){
  double cdf_normal = 0.0;
  
  cdf_normal = cdf(normal_mdist(mean,sd), q);
  
  if (cdf_normal == 0) {
    double dbl_m = std::numeric_limits< double >::min();
    cdf_normal = dbl_m;
  }
  
  return(cdf_normal);
}

/* Exponential PDF + CDF */
/* MEAN = 1/RATE */
inline double pdf_exponential(double x, double rate) {
  double pdf_exp = 0.0;
  pdf_exp = pdf(exp_mdist(rate), x);
  
  if (pdf_exp == 0) {
    double dbl_min = std::numeric_limits< double >::min();
    pdf_exp = dbl_min;
  }
  
  return(pdf_exp);
}

inline double cdf_exponential(double q, double rate) {
  double cdf_exp = 0.0;
  cdf_exp = cdf(exp_mdist(rate), q);
  
  if (cdf_exp == 0) {
    double dbl_m = std::numeric_limits< double >::min();
    cdf_exp = dbl_m;
  }
  
  return(cdf_exp);
}

#endif