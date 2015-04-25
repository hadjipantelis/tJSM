
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]


Eigen::MatrixXd calc_M1_M2_Hadamard_y2( Eigen::Map<Eigen::ArrayXXd> & A1, const Eigen::Map<Eigen::ArrayXXd> & A2, const Eigen::Map<Eigen::VectorXd> & v3, const int u){ 
  
  //  This function implements:
  //  do.call(cbind, lapply(1:u, function(i) as.vector((A1 * A2[,i]) %*% v3))) # N*u matrix #
  
  unsigned int N = A1.cols();

  Eigen::MatrixXd Mr = Eigen::MatrixXd( A1.rows(), u);  
  for(int j=0; j != u; ++j){
   
  Eigen::ArrayXXd Ak = A1; 
    for(unsigned int i=0; i != N; ++i){
      Ak.col(i)  *= A2.col(j) ;
    } 
    Mr.col(j)  = Ak.matrix() * v3 ;
  } 
  return (  Mr  );  
}
 
