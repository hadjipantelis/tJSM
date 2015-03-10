example_yT_Minv_y <- cxxfunction( plugin = "RcppEigen",  signature(y_i = "vector", M_i ="matrix"), body='	
	#include <iostream> 
 	using Eigen::Map;		 			// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 						// to do the LLT decomposition 


	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	//Map vector M_i to vectorXd theta /Hyperparameters Theta
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	//Map vector y_i to vectorXd y / Measurements Values 
 	LLT<MatrixXd> llt_M(M);
	//VectorXd alpha = llt_M.solve(y);  
	return Rcpp::wrap( llt_M.solve(y).transpose() * y );  
	 
' )

example_yT_Minv  <- cxxfunction( plugin = "RcppEigen",  signature(y_i = "vector", M_i ="matrix"), body='	
	#include <iostream> 
 	using Eigen::Map;		 			// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 						// to do the LLT decomposition 


	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	//Map vector M_i to vectorXd theta /Hyperparameters Theta
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	//Map vector y_i to vectorXd y / Measurements Values 
 	LLT<MatrixXd> llt_M(M);
	//VectorXd alpha = llt_M.solve(y);  
	return Rcpp::wrap( llt_M.solve(y) );  
	 
' )


%A.triangularView<Upper>().solve<OnTheRight>(bT)

example_inv_cholM <- cxxfunction( plugin = "RcppEigen",  signature(M_i ="matrix"), body='	
	#include <iostream> 
 	using Eigen::Map;		 			// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::LLT; 						// to do the LLT decomposition  

	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	//Map vector M_i to vectorXd theta /Hyperparameters Theta

	int p = M.cols();
	MatrixXd b = MatrixXd::Identity(p, p);
 	//std::cout << MatrixXd( m.llt().solve(b).llt().matrixU()) << std::endl;
	return Rcpp::wrap( MatrixXd( M.llt().solve(b).llt().matrixU()) );  
	 
' )

example_expM <- cxxfunction( plugin = "RcppEigen",  signature(M_i ="array"), body='	
	#include <iostream> 
 	using Eigen::Map;		 			// to map input variable to an existing array of data 
	using Eigen::ArrayXd; 				// to use MatrixXd 

	const Map<ArrayXd> M(Rcpp::as<Map<ArrayXd> > (M_i));	//Map vector M_i to vectorXd theta /Hyperparameters Theta

	//std::cout << MatrixXd( m.llt().solve(b).llt().matrixU()) << std::endl;
	return Rcpp::wrap( (M.exp()) );   
' )


example_Mymult <- cxxfunction( plugin = "RcppEigen", signature(y_i = "vector", M_i ="matrix"), body='	
	#include <iostream> 
 	using Eigen::Map;		 			// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd

	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	//Map vector M_i to vectorXd theta /Hyperparameters Theta
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	//Map vector y_i to vectorXd y / Measurements Values 

	return Rcpp::wrap(  M * y  );  
	 
' )


example_MMMmult <- cxxfunction( plugin = "RcppEigen",  signature(M_i1 ="matrix", M_i2 ="matrix", M_i3 ="matrix"), body='	
	#include <iostream> 
 	using Eigen::Map;		 			// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd 

	const Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));	//Map vector M_i to vectorXd theta /Hyperparameters Theta
	const Map<MatrixXd> M2(Rcpp::as<Map<MatrixXd> > (M_i2));	//Map vector M_i to vectorXd theta /Hyperparameters Theta
	const Map<MatrixXd> M3(Rcpp::as<Map<MatrixXd> > (M_i3));	//Map vector M_i to vectorXd theta /Hyperparameters Theta

	//std::cout << MatrixXd( m.llt().solve(b).llt().matrixU()) << std::endl;
	return Rcpp::wrap( M1.cwiseProduct(M2).cwiseProduct(M3)  );  
	 
' )

 
 



#include <Rcpp.h>
using namespace Rcpp;

#include <cstdlib>
using std::size_t;

void hadamardMultiplyMatrixByVectorInPlace(double* restrict x,
                                           size_t numRows, size_t numCols,
                                           const double* restrict y)
{
  if (numRows == 0 || numCols == 0) return;

  for (size_t col = 0; col < numCols; ++col) {
    double* restrict x_col = x + col * numRows;

    for (size_t row = 0; row < numRows; ++row) {
      x_col[row] *= y[row];
    }
  }
}

// [[Rcpp::export]]
NumericMatrix C_matvecprod_elwise_inplace(NumericMatrix& X,
                                          const NumericVector& y)
{
  // do some dimension checking here

  hadamardMultiplyMatrixByVectorInPlace(X.begin(), X.nrow(), X.ncol(),
                                        y.begin());

  return X;
}















