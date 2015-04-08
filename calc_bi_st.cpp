
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::MatrixXd calc_MVND(const Eigen::Map<Eigen::VectorXd> & x, const Eigen::Map<Eigen::VectorXd> & mu, const Eigen::Map<Eigen::MatrixXd> & K){ 
 
  Eigen::LLT<Eigen::MatrixXd> llt_M(M1.inverse()); 
  Eigen::MatrixXd Res = sqrt(2.) *llt_M.matrixLLT().triangularView<Eigen::Lower>().transpose().solve(M2.transpose()); 
 
  for (unsigned int i = 0; i != Res.cols(); ++i) {
    Res.col(i) += y0;    
  }
return Res;



cat('Compliling calc_MVND\n')
calc_MVND <- cxxfunction( plugin = "RcppEigen",  signature( y_i1 = "vector", y_i2 = "vector", M_i ="matrix"), body='
  // Calculate $MVND IN 2D$
  using Eigen::Map;       // to map input variable to an existing array of data
  using Eigen::MatrixXd;
  using Eigen::VectorXd;       // to use VectorXd
  using Eigen::LLT;       // to do the LLT decomposition  
  using Eigen::ArrayXd;
  using Eigen::Lower;   
  const Map<MatrixXd> K(Rcpp::as<Map<MatrixXd> > (M_i));  // Map matrix M_i to matrixXd M
  const Map<VectorXd> x(Rcpp::as<Map<VectorXd> > (y_i1));  // Map vector y_i to vectorXd y  
  const Map<VectorXd> mu(Rcpp::as<Map<VectorXd> > (y_i2));// Map vector y_i to vectorXd y  
 
  double p = x.size();
  LLT<MatrixXd> LLT_of_K(K); // compute the Cholesky decomposition of K
  if ( !LLT_of_K.info() ) { 
    MatrixXd Rooti = LLT_of_K.matrixLLT().triangularView<Lower>().solve(MatrixXd::Identity(p,p)); 
    double quads = (Rooti * (x-mu)).array().square().sum();
    double res   =  exp( -((p/2.0)*log(2.0*M_PI))  + Rooti.diagonal().array().log().sum() - 0.5*quads); 
    return Rcpp::wrap( res );
  } else {
    return Rcpp::wrap( "Cholesky decomposition failed" );
  }
' )
} 
