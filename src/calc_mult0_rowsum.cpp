
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
Eigen::MatrixXd calc_mult0_rowsum(const Eigen::Map<Eigen::VectorXi> & v, const Eigen::Map<Eigen::VectorXd> & u, const Eigen::Map<Eigen::MatrixXd> & M){ 
         
  const unsigned int l = v.size();
  const unsigned int m = M.cols();  

  Eigen::MatrixXd Res = Eigen::MatrixXd::Zero(v.maxCoeff() ,m);  

  for (unsigned int i = 0; i < m; ++i){
  unsigned int k = 0;
    for(unsigned int j = 0; j != l; ++j){
     Res(k,i) = Res(k,i) + M(j,i) * u(j); 
      if (j  < M.rows() ) {
        if( v[j] != v[j+1] ){ 
          k++;
        } 
      }
    }
  }

  return( Res );    
}

