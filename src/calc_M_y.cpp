
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::VectorXd calc_M_y( const Eigen::Map<Eigen::VectorXd> & v, const Eigen::Map<Eigen::MatrixXd> & M){ 
  return(M * v);  
}


