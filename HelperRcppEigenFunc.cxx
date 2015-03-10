library(inline,RcppEigen)

# Modify the plugin for RcppEigen to support OpenMP
settingsE <- getPlugin("RcppEigen")
#//settingsE$env$PKG_CXXFLAGS <- paste('-fopenmp', settingsE$env$PKG_CXXFLAGS)
#//settingsE$env$PKG_LIBS <- paste('-fopenmp -lgomp', settingsE$env$PKG_LIBS)


print('Compliling calc_yT_Minv_y')
calc_yT_Minv_y <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(y_i = "vector", M_i ="matrix"), body='
	// Calculate $y^T M^{-1} y$	
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition 

	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	// Map matrix M_i to matrixXd M  
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	// Map vector y_i to vectorXd y  
 	LLT<MatrixXd> llt_M(M);					// Calculate the LLT decomposition 

	return Rcpp::wrap( llt_M.solve(y).transpose() * y );  	 
' )

print('Compliling calc_yT_Minv')
calc_yT_Minv  <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(y_i = "vector", M_i ="matrix"), body='	
	// Calculate $y^T M^{-1} $	
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition 

	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	// Map matrix M_i to matrixXd M  
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	// Map vector y_i to vectorXd y  
 	LLT<MatrixXd> llt_M(M); 				// Calculate the LLT decomposition 

	return Rcpp::wrap( llt_M.solve(y) );  	 
' )

print('Compliling calc_chol_Minv')
calc_chol_Minv <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i ="matrix"), body='	
	// Calculate $chol(M^{-1})$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::LLT; 				// to do the LLT decomposition  

	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	// Map matrix M_i to matrixXd M
	unsigned int p = M.cols(); 				// Get the number of columns in M
	MatrixXd B = MatrixXd::Identity(p, p); 			// Define an identity matrix B

	return Rcpp::wrap( MatrixXd( M.llt().solve(B).llt().matrixU()) );  
	 
' )

print('Compliling calc_expM')
calc_expM <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i ="array"), body='	
	// Calculate $exp(M)$ - this function returns an array
 	using Eigen::Map;		 		// to map input variable to an existing array of data 
	using Eigen::ArrayXd; 				// to use ArrayXd 

	const Map<ArrayXd> A(Rcpp::as<Map<ArrayXd> > (M_i));	//Map vector M_i to arrayXd A  
 
	return Rcpp::wrap( (A.exp()) );   
' )


print('Compliling calc_M_y')
calc_M_y <- cxxfunction(settings=settingsE, plugin = "RcppEigen", signature(y_i = "vector", M_i ="matrix"), body='	
	// Calculate $M y$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd

	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	//Map vector M_i to MatrixXd M 
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	//Map vector y_i to vectorXd y 

	return Rcpp::wrap(  M * y  );  	 
' )


print('Compliling calc_M1_M2_M3')
calc_M1_M2_M3 <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 ="matrix", M_i2 ="matrix", M_i3 ="matrix"), body='	
	// Calculate the Hadamard product $M1 M2 M3$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd 

	const Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));	//Map vector M_1 to MatrixXd M1
	const Map<MatrixXd> M2(Rcpp::as<Map<MatrixXd> > (M_i2));	//Map vector M_2 to MatrixXd M2
	const Map<MatrixXd> M3(Rcpp::as<Map<MatrixXd> > (M_i3));	//Map vector M_3 to MatrixXd M3
 
	return Rcpp::wrap( M1.cwiseProduct(M2).cwiseProduct(M3)  );  	 
' )

print('Compliling calc_M1_M2_M3_Hadamard')
calc_M1_M2_M3_Hadamard <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 ="matrix", M_i2 ="matrix", M_i3 ="matrix", v_i = "vector"), body='	
	// Calculate the Hadamard product $M1 M2 M3$ using indeces at v
	// This function does IN-PLACE multiplication
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd 
	using Eigen::VectorXi; 				// to use VectorXd 

	Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));	//Map matrix M_1 to MatrixXd M1
	const Map<MatrixXd> M2(Rcpp::as<Map<MatrixXd> > (M_i2));//Map matrix M_2 to MatrixXd M2
	const Map<MatrixXd> M3(Rcpp::as<Map<MatrixXd> > (M_i3));//Map matrix M_3 to MatrixXd M3
	const Map<VectorXi> v(Rcpp::as<Map<VectorXi> > (v_i));	//Map vector v_i to MatrixXd v
	unsigned int N = v.size();

	for(unsigned int i=0; i != N; ++i){
		 M1.row(i).array()  *= M2.row(v(i)).array() *  M3.row(v(i)).array();
	}
 
	//return Rcpp::wrap( M1 ); // This function returns VOID
' )

print('Compliling calc_y_M1_Hadamard')
calc_y_M_Hadamard <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 ="matrix", v_i1 = "vector", v_i2 = "vector"), body='	
	// Calculate the product $ y$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd 
	using Eigen::VectorXd; 				// to use VectorXd 
	using Eigen::VectorXi; 				// to use VectorXd 

	Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));		//Map matrix M_1 to MatrixXd M1
	const Map<VectorXd> v1(Rcpp::as<Map<VectorXd> > (v_i1));	//Map vector v_i to MatrixXd v
	const Map<VectorXi> v2(Rcpp::as<Map<VectorXi> > (v_i2));	//Map vector v_i to MatrixXd v
	unsigned int N = v2.size();

	for(unsigned int i=0; i != N; ++i){
		 M1.row(i).array()  *=  v1(v2(i));
	}	 
' )


 
