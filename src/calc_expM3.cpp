
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
void calc_expM3(Eigen::Map<Eigen::ArrayXXd> & A){ 
  A = (A.exp()); 
}

