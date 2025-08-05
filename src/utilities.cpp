#include <Rcpp.h>
#include "utilities.h"
using namespace Rcpp;

// [[Rcpp::export]]
double mu(double t, double wp, double tp, double dp, double wr){
  if(t > (tp - wp) && t <= tp){
    return((dp/wp)*(t-(tp-wp)));
  } else if(t > tp && t <= (tp + wr)){
    return(dp - (dp/wr)*(t - tp));
  } else{
    return(0.0);
  }
}