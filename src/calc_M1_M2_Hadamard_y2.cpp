
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]


Eigen::MatrixXd calc_M1_M2_Hadamard_y2( Eigen::Map<Eigen::ArrayXXd> & A1, const Eigen::Map<Eigen::ArrayXXd> & A2, const Eigen::Map<Eigen::VectorXd> & v3, const int u){ 
  
  // Calculate the Hadamard product $M_i1 M_i2 M_i3$ using indeces at v_i
  // This function makes in-place computations

  unsigned int N = A1.cols();

  Eigen::MatrixXd Mr = Eigen::MatrixXd( A1.rows(), A1.cols());  
  for(int j=0; j != u; ++j){
  Eigen::ArrayXXd Ak = A1; 
    for(unsigned int i=0; i != N; ++i){
      Ak.col(i)  *= A2.col(j) ;
    } 
     Mr.col(j)  = Ak.matrix() * v3 ;
  } 
  return (  Mr  );  
}
 
