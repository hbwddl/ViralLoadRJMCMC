#include <Rcpp.h>
#include "debug.h"
using namespace Rcpp;

void print_pos(std::string file_arg,
               int line_arg){
  Rcout << file_arg << " line " << line_arg << "\n";
}
