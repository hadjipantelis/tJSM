
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
void calc_expM2(Eigen::Map<Eigen::ArrayXd> & A){ 
  A = (A.exp()); 
}

