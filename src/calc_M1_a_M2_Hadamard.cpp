
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]


void calc_M1_a_M2_Hadamard(Eigen::Map<Eigen::MatrixXd> & M1, const Eigen::Map<Eigen::MatrixXd> & M2, const double a, const Eigen::Map<Eigen::VectorXi> & v){ 
  
  // Calculate the Hadamard product $M_i1 M_i2 M_i3$ using indeces at v_i
  // This function makes in-place computations

  unsigned int N = v.size();
  
  M1.array() *= a;

  for(unsigned int i=0; i != N; ++i){
    M1.row(i).array()  *= M2.row(v(i)).array();
  } 

}
