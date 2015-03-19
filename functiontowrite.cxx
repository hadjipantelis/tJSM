
if(11==1){
print('Compliling calc_bi_st')
calc_bi_st <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature( y_i0 = "numeric", y_i1 = "matrix", M_i ="matrix"), body='
	// Calculate $chol(M_i^{-1})$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition  
	using Eigen::Lower;
	using Eigen::Upper;
#include <iostream>
 
	const Map<VectorXd> y0(Rcpp::as<Map<VectorXd> > (y_i0));	// Map vector y_i to vectorXd y  
	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));		// Map matrix M_i to matrixXd M
	const Map<MatrixXd> y1(Rcpp::as<Map<MatrixXd> > (y_i1));	// Map vector y_i to vectorXd y  
	
	unsigned int p = M.cols();
	LLT<MatrixXd> llt_M(M.inverse());  
return Rcpp::wrap( sqrt(2.) * llt_M.matrixLLT().triangularView<Lower>().transpose().solve(y1));
' )


}
if(1==1){
print('Compliling calc_bi_st')
calc_solveEig <- cxxfunction( plugin = "RcppEigen",  signature(y_i0 = "vector", y_i1 = "matrix", M_i ="matrix"), body='
	// Calculate $chol(M_i^{-1})$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition  
	using Eigen::Lower;
	using Eigen::Upper;
	#include <iostream>
  
	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));		// Map matrix M_i to matrixXd M
	const Map<MatrixXd> y1(Rcpp::as<Map<MatrixXd> > (y_i1));	// Map vector y_i to vectorXd y  
	const Map<VectorXd> y0(Rcpp::as<Map<VectorXd> > (y_i0));	// Map vector y_i to vectorXd y  
	
	unsigned int p = M.cols();
	LLT<MatrixXd> llt_M(M.inverse()); 
	MatrixXd Res = sqrt(2.) *llt_M.matrixLLT().triangularView<Lower>().transpose().solve(y1.transpose());
	//std::cout << Res << std::endl;
	//std::cout << Res.rows() << std::endl;
	//std::cout << Res.cols() << std::endl;
 	for (unsigned int i = 0; i != Res.cols(); ++i) {
		Res.col(i) += y0;		
	}

return Rcpp::wrap( Res  );
' )
}
