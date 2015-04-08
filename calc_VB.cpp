#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::MatrixXd calc_VB( const Eigen::Map<Eigen::MatrixXd> & M1, const Eigen::Map<Eigen::MatrixXd> & M2, const Eigen::Map<Eigen::MatrixXd> & M3){ 

  Eigen::LLT<Eigen::MatrixXd> llt_M3(M3);          // Calculate the LLT decomposition 
  
  Eigen::MatrixXd b = M2.transpose() * llt_M3.solve(M2) ; 
  
  return( M1 - M1 * b * M1 );  
}



