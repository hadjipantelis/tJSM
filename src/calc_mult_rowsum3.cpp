
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
Rcpp::List calc_mult_rowsum3(const Eigen::Map<Eigen::ArrayXi> & v, const Eigen::Map<Eigen::ArrayXXd> & B, const Eigen::Map<Eigen::ArrayXXd> & M, const Eigen::Map<Eigen::ArrayXd> & A, const double ncb2){ 
  
  const unsigned int n =  uint(ncb2);
  const unsigned int l = v.size();
  const unsigned int m = M.cols();
       
  Rcpp::List output(n);

 for (unsigned int u = 0; u != n; u++){
    Eigen::ArrayXXd Res = Eigen::ArrayXXd::Zero(v.maxCoeff() ,m);  
    for (unsigned int i = 0; i < m; ++i){
      unsigned int k = 0;
      for(unsigned int j = 0; j != l; ++j){
        Res(k,i) = Res(k,i) + M(j,i) * B(j,u); 
        if (j  < M.rows() ) {
          if( v[j] != v[j+1] ){ 
            k++;
          } 
        }
      }
    }
    Res.array() *= A;
    output[u] = Res;
  }
  return( output );    
}

