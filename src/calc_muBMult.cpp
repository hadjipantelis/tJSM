
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::VectorXd calc_muBMult( const Eigen::Map<Eigen::MatrixXd> & BSold, const Eigen::Map<Eigen::MatrixXd> & VY, const Eigen::Map<Eigen::VectorXd> & BTg, const Eigen::Map<Eigen::VectorXd> & Yst){ 

  Eigen::VectorXd yf = Yst - BTg;      
  Eigen::LLT<Eigen::MatrixXd> llt_VY(VY);         // Calculate the LLT decomposition 
        
  return (BSold * BTg.transpose() * llt_VY.solve(yf) );    
 
}
