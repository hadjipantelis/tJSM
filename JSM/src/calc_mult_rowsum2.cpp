#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
Eigen::MatrixXd calc_mult_rowsum2(const Eigen::Map<Eigen::VectorXi> & v, const Eigen::Map<Eigen::MatrixXd> & L, const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::ArrayXd> & A){ 

  //  This function implements:
  //  A * rowsum( L * M  , v))
  //  as before 'v' needs to be sorted!!
        
  const unsigned int l = v.size();
  const unsigned int mc = M.cols();  
  const unsigned int mr = M.rows();

  Eigen::MatrixXd Res = Eigen::MatrixXd::Zero(v.maxCoeff() ,mc);  

  for (unsigned int i = 0; i < mc; ++i){
  unsigned int k = 0;
    for(unsigned int j = 0; j != l; ++j){
     Res(k,i) = Res(k,i) + M(j,i) * L(j,i); 
      if (j  < mr ) {
        if( v[j] != v[j+1] ){ 
          k++;
        } 
      }
    }
  }

  Res.array() *= A;

  return( Res );    
}
