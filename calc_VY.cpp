#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::MatrixXd calc_VY( const Eigen::Map<Eigen::MatrixXd> & v, const Eigen::Map<Eigen::MatrixXd> & a, const double b){ 

 Eigen::MatrixXd M1 = v*a*v.transpose();
  M1.diagonal().array() += b;
  return M1;
}


