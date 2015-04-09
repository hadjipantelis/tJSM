
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]


void calc_y_a(Eigen::Map<Eigen::ArrayXd> & v, const double & a){

  // Calculate the vector multiplication $a v_i1$ a_i being a scalar
  // This function makes in-place computations 

  v *= a;  

}
