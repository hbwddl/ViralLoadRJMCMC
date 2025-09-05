#include <Rcpp.h>
#include "debug.h"
using namespace Rcpp;

void print_pos(std::string file_arg,
               int line_arg,
               int debug_true){
  if(debug_true == 1){
    Rcout << file_arg << " line " << line_arg << "\n";
  }
}
