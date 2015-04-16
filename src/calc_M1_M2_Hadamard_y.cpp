
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]


Eigen::VectorXd calc_M1_M2_Hadamard_y(Eigen::Map<Eigen::MatrixXd> & M1, const Eigen::Map<Eigen::ArrayXd> & v1, const Eigen::Map<Eigen::VectorXd> & v2){ 
  
  // Calculate the Hadamard product $M_i1 M_i2 M_i3$ using indeces at v_i
  // This function makes in-place computations

  unsigned int N = M1.cols();

  Eigen::MatrixXd Mk = M1; 

  for(unsigned int i=0; i != N; ++i){
     Mk.col(i).array()  *= v1 ;
  } 
  
  return (  Mk * v2  );  

}
