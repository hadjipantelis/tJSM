library(inline,RcppEigen)

# Modify the plugin for RcppEigen to support OpenMP
settingsE <- getPlugin("RcppEigen")
#//settingsE$env$PKG_CXXFLAGS <- paste('-fopenmp', settingsE$env$PKG_CXXFLAGS)
#//settingsE$env$PKG_LIBS <- paste('-fopenmp -lgomp', settingsE$env$PKG_LIBS)


print('Compliling calc_yT_Minv_y')
calc_yT_Minv_y <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(y_i = "vector", M_i ="matrix"), body='
	// Calculate $y_i^T M_i^{-1} y$	
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
	// Calculate $y_i^T M_i^{-1} $	
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
	// Calculate $chol(M_i^{-1})$
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
	// Calculate $exp(M_i)$ - this function returns an array
 	using Eigen::Map;		 		// to map input variable to an existing array of data 
	using Eigen::ArrayXd; 				// to use ArrayXd 

	const Map<ArrayXd> A(Rcpp::as<Map<ArrayXd> > (M_i));	//Map vector M_i to arrayXd A  
 
	return Rcpp::wrap( (A.exp()) );   
' )


print('Compliling calc_M_y')
calc_M_y <- cxxfunction(settings=settingsE, plugin = "RcppEigen", signature(y_i = "vector", M_i ="matrix"), body='	
	// Calculate $M_i y_i$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd

	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	//Map vector M_i to MatrixXd M 
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	//Map vector y_i to vectorXd y 

	return Rcpp::wrap(  M * y  );  	 
' )

print('Compliling calc_M1_M2_M3')
calc_M1_M2_M3 <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 ="matrix", M_i2 ="matrix", M_i3 ="matrix"), body='	
	// Calculate the Hadamard product $M_i1 M_i2 M_i3$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd 

	const Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));	//Map vector M_i1 to MatrixXd M1
	const Map<MatrixXd> M2(Rcpp::as<Map<MatrixXd> > (M_i2));	//Map vector M_i2 to MatrixXd M2
	const Map<MatrixXd> M3(Rcpp::as<Map<MatrixXd> > (M_i3));	//Map vector M_i3 to MatrixXd M3
 
	return Rcpp::wrap( M1.cwiseProduct(M2).cwiseProduct(M3)  );  	 
' )

print('Compliling calc_M1_M2_M3_Hadamard')
calc_M1_M2_M3_Hadamard <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 ="matrix", M_i2 ="matrix", M_i3 ="matrix", v_i = "vector"), body='	
	// Calculate the Hadamard product $M_i1 M_i2 M_i3$ using indeces at v_i
	// This function does IN-PLACE multiplication
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd 
	using Eigen::VectorXi; 				// to use VectorXi 

	Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));	//Map matrix M_1 to MatrixXd M1
	const Map<MatrixXd> M2(Rcpp::as<Map<MatrixXd> > (M_i2));//Map matrix M_2 to MatrixXd M2
	const Map<MatrixXd> M3(Rcpp::as<Map<MatrixXd> > (M_i3));//Map matrix M_3 to MatrixXd M3
	const Map<VectorXi> v(Rcpp::as<Map<VectorXi> > (v_i));	//Map vector v_i to VectorXi v (integer)
	unsigned int N = v.size();

	for(unsigned int i=0; i != N; ++i){
		 M1.row(i).array()  *= M2.row(v(i)).array() *  M3.row(v(i)).array();
	}
 
	//return Rcpp::wrap( M1 ); // This function returns VOID
' )

print('Compliling calc_y_M1_Hadamard')
calc_y_M_Hadamard <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 ="matrix", v_i1 = "vector", v_i2 = "vector"), body='	
	// Calculate the Hadamard product $ M_i1 v_i1$ using indeces at v_i2
	// This function does IN-PLACE multiplication
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd 
	using Eigen::VectorXd; 				// to use VectorXd 
	using Eigen::VectorXi; 				// to use VectorXi

	Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));		//Map matrix M_1 to MatrixXd M1
	const Map<VectorXd> v1(Rcpp::as<Map<VectorXd> > (v_i1));	//Map vector v_i to VectorXd v1
	const Map<VectorXi> v2(Rcpp::as<Map<VectorXi> > (v_i2));	//Map vector v_i to VectorXi v2
	unsigned int N = v2.size();

	for(unsigned int i=0; i != N; ++i){
		 M1.row(i).array()  *=  v1(v2(i));
	}
 
	//return Rcpp::wrap( M1 ); // This function returns VOID	 
' )

print('Compliling calc_y_a')
calc_y_a <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(v_i1 = "array", a_i = "numeric"), body='	
	// Calculate the vector multiplication $a v_i1$ a_i being a scalar
	// This function does IN-PLACE multiplication
 	using Eigen::Map;		 		// to map input variable to an existing array of data 
	using Eigen::ArrayXd; 				// to use ArrayXd 

	Map<ArrayXd> v1(Rcpp::as<Map<ArrayXd> > (v_i1));	//Map array v_i to ArrayXd v 
	const double a (Rcpp::as<double> (a_i));		//Map double a_i to double a 
 
	v1 *= a; 

	//return Rcpp::wrap( v1 ); // This function returns VOID	
' )


print('Compliling calc_tapply_vect_sum')
calc_tapply_vect_sum <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(v_i1 = "array", v_i2 = "array"), body='	
	// Calculate the equivalent of tapply(v_i1, v_i2, sum)
 	using Eigen::Map;		 		// to map input variable to an existing array of data 
	using Eigen::ArrayXd; 				// to use ArrayXd 
	using Eigen::ArrayXi; 				// to use ArrayXi (integers)  

	const 	Map<ArrayXd> v1(Rcpp::as<Map<ArrayXd> > (v_i1));		//Map array v_i1 to ArrayXd v1 
	const 	Map<ArrayXi> v2(Rcpp::as<Map<ArrayXi> > (v_i2));		//Map array v_i2 to ArrayXi v2 
	unsigned int N = v2.size();
	unsigned int a = v2.maxCoeff();
	ArrayXd v3 = ArrayXd::Zero(a+1);	
	
	for( unsigned int i=0; i !=N; i++){	
		v3(v2(i)) += v1(i);		
	} 

	return Rcpp::wrap( v3 );  
' )


print('Compliling calc_VY')
calc_VY <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(v_i1 = "vector", a_i = "numeric", b_i = "numeric"), body='	
	// Calculate VY as $v_i * a_i * v_i^T + I * b_i$ 
 	using Eigen::Map;		 		// to map input variable to an existing array of data 
	using Eigen::ArrayXd; 				// to use ArrayXd  
	using Eigen::MatrixXd; 				// to use MatrixXd 
	using Eigen::VectorXd; 				// to use VectorXd

	const 	Map<VectorXd> v(Rcpp::as<Map<VectorXd> > (v_i1));	//Map array v_i1 to ArrayXd v1  
	const double a (Rcpp::as<double> (a_i));			//Map double a_i to double a 
	const double b (Rcpp::as<double> (b_i));			//Map double b_i to double b 

	MatrixXd M1 = v*a*v.transpose();
	M1.diagonal().array() += b;

	return Rcpp::wrap( M1 );  
' )


print('Compliling calc_VB')
calc_VB <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(y_i = "vector", M_i ="matrix", a_i = "numeric"), body='
	// Calculate VB  as   a_i  -  a_i * y_i^T * M_i^(-1) * y_i * a_i
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition 

	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	// Map matrix M_i to matrixXd M  
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	// Map vector y_i to vectorXd y  
	const double a (Rcpp::as<double> (a_i));		// Map double a_i to double a 
 	LLT<MatrixXd> llt_M(M);					// Calculate the LLT decomposition 
	
	double b = a; 
	b = (llt_M.solve(y).transpose() * y); 
	return Rcpp::wrap( a - b * a * a   );  	 
' )


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


print('Compliling calc_bi_st')
calc_bi_st <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature( a_i = "numeric", y_i = "vector", M_i ="matrix"), body='
	// Calculate $chol(M_i^{-1})$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition  
	using Eigen::Lower;

	const double a (Rcpp::as<double> (a_i));		// Map double a_i to double a 
	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	// Map matrix M_i to matrixXd M
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	// Map vector y_i to vectorXd y  
	// unsigned int p = M.cols(); 				// Get the number of columns in M
	// MatrixXd B = MatrixXd::Identity(p, p); 			// Define an identity matrix B 
  	return Rcpp::wrap( a +  sqrt ( 2.*  M(0,0) ) * y.transpose().array());
	//return Rcpp::wrap( (  1.414213562373095*M.llt().solve(B).llt().matrixLLT().triangularView<Lower>().solve(y)) )
' )

print('Compliling calc_MVND')
calc_MVND <- cxxfunction( plugin = "RcppEigen",  signature( y_i1 = "vector", y_i2 = "vector", M_i ="matrix"), body='
	// Calculate $MVND IN 2D$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd;
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition  
	using Eigen::ArrayXd;
	using Eigen::Lower;	 
	const Map<MatrixXd> K(Rcpp::as<Map<MatrixXd> > (M_i));	// Map matrix M_i to matrixXd M
	const Map<VectorXd> x(Rcpp::as<Map<VectorXd> > (y_i1));	// Map vector y_i to vectorXd y  
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

print('Compliling calc_rowsum')
calc_rowsum <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(y_i = "vector", M_i = "Matrix"), body='
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXi; 				// to use VectorXd
	//#include <iostream>

	const Map<VectorXi> y1(Rcpp::as<Map<VectorXi> > (y_i));	// Map vector y_i1 to vectorXd y1  
	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));		// Map matrix M_i to vectorXd M  
	
	const unsigned int l = y1.size();
	const unsigned int m = M.cols(); 	

	MatrixXd Res = MatrixXd::Zero(y1.maxCoeff() ,m);	

	for (unsigned int i = 0; i < m; ++i){
	unsigned int k = 0;
		for(unsigned int j = 0; j != l; ++j){
 			Res(k,i) = Res(k,i) + M(j,i); 
			if (j  < M.rows() ) {
				if( y1[j] != y1[j+1] ){ 
					k++;
				} 
			}
		}
	}
	return Rcpp::wrap( Res );  	 
' )

