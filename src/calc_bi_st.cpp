
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::MatrixXd calc_bi_st (const Eigen::Map<Eigen::VectorXd> & v0, const Eigen::Map<Eigen::VectorXd> & v1, const Eigen::Map<Eigen::MatrixXd> & M){ 
 
  Eigen::LLT<Eigen::MatrixXd> llt_M(M.inverse()); 
  Eigen::MatrixXd Res = sqrt(2.) *llt_M.matrixLLT().triangularView<Eigen::Lower>().transpose().solve(v1.transpose()); 
 
  for (unsigned int i = 0; i != Res.cols(); ++i) {
    Res.col(i) += v0;    
  }
return Res;
 
}  
