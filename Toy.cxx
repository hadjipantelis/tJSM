
print('Compliling calc_VB')
calc_VB <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(y_i = "vector", M_i ="matrix", a_i = "numeric"), body='
	// Calculate $y_i^T M_i^{-1} y$	
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition 

	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	// Map matrix M_i to matrixXd M  
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	// Map vector y_i to vectorXd y  
	const double a (Rcpp::as<double> (a_i));			//Map vector v_i to MatrixXd v 
 	LLT<MatrixXd> llt_M(M);					// Calculate the LLT decomposition 
	
	double b = a; 
	b = (llt_M.solve(y).transpose() * y); 
	return Rcpp::wrap( a - b * a* a   );  	 
' )

