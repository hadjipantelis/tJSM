print('Compliling calc_muB')
calc_muB <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature( a_i = "numeric", y_i0 = "vector", y_i1 = "vector", y_i2 = "vector", M_i1 ="matrix", M_i2 ="matrix"), body='
	//  as.vector(BSigma.old %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - X.st[[i]] %*% beta.old)))
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition  
	using Eigen::Lower;

	const double bsigmaold (Rcpp::as<double> (a_i));		// Map double a_i to double a 
	const Map<MatrixXd> VY(Rcpp::as<Map<MatrixXd> > (M_i1));	// Map matrix M_i1 to matrixXd M
	const Map<VectorXd> Zst(Rcpp::as<Map<VectorXd> > (y_i0));	// Map vector y_i0 to vectorXd y  
	const Map<MatrixXd> Xst(Rcpp::as<Map<MatrixXd> > (M_i2));	// Map matrix M_i2 to matrixXd M
	const Map<VectorXd> Yst(Rcpp::as<Map<VectorXd> > (y_i1));	// Map vector y_i1 to vectorXd y  
	const Map<VectorXd> betaold(Rcpp::as<Map<VectorXd> > (y_i2));	// Map vector y_i3 to vectorXd y  

	VectorXd yf = Yst - Xst * betaold;			
 	LLT<MatrixXd> llt_VY(VY); 					// Calculate the LLT decomposition 

	return Rcpp::wrap( bsigmaold * Zst.transpose() * llt_VY.solve(yf) );  	 
' )



# ASDF <- calc_muB( BSigma.old,  y_i0=Z.st[[i]], y_i1=Y.st[[i]],  y_i2=beta.old, M_i1=VY[[i]], M_i2=X.st[[i]])
