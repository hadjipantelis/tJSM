
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
void calc_expM3(Eigen::Map<Eigen::ArrayXXd> & A){ 
  
  //  This function implements:
  //  A = exp(A);
  //  This function makes in-place computations

  A = (A.exp()); 
}

