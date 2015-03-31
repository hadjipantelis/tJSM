library(inline,RcppEigen)

# Modify the plugin for RcppEigen to support OpenMP
settingsE <- getPlugin("RcppEigen")
#//settingsE$env$PKG_CXXFLAGS <- paste('-fopenmp', settingsE$env$PKG_CXXFLAGS)
#//settingsE$env$PKG_LIBS <- paste('-fopenmp -lgomp', settingsE$env$PKG_LIBS)

# These functions were used mostly for testing or middlemen
if(1==2){ 
cat('Compliling calc_yT_Minv_y\n')
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

cat('Compliling calc_yT_Minv\n')
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

cat('Compliling calc_chol_Minv\n')
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

cat('Compliling calc_y_M1_Hadamard\n')
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


}

cat('Compliling calc_expM\n')
calc_expM <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i ="array"), body='	
	// Calculate $exp(M_i)$ - this function returns an array
 	using Eigen::Map;		 		// to map input variable to an existing array of data 
	using Eigen::ArrayXd; 				// to use ArrayXd 

	const Map<ArrayXd> A(Rcpp::as<Map<ArrayXd> > (M_i));	//Map vector M_i to arrayXd A  
 
	return Rcpp::wrap( (A.exp()) );   
' )


cat('Compliling calc_M_y\n')
calc_M_y <- cxxfunction(settings=settingsE, plugin = "RcppEigen", signature(y_i = "vector", M_i ="matrix"), body='	
	// Calculate $M_i y_i$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd

	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));	//Map vector M_i to MatrixXd M 
	const Map<VectorXd> y(Rcpp::as<Map<VectorXd> > (y_i));	//Map vector y_i to vectorXd y 

	return Rcpp::wrap(  M * y  );  	 
' )

cat('Compliling calc_M1_M2_Hadamard\n')
calc_M1_M2_Hadamard <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 ="array", M_i2 ="array" ), body='	
	// Calculate the Hadamard product $M_i1 M_i2 M_i3$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::ArrayXd; 				// to use MatrixXd 

	Map<ArrayXd> M1(Rcpp::as<Map<ArrayXd> > (M_i1));	//Map vector M_i1 to MatrixXd M1
	const Map<ArrayXd> M2(Rcpp::as<Map<ArrayXd> > (M_i2));	//Map vector M_i2 to MatrixXd M2 
  
		 M1 *= M2; 

	// return Rcpp::wrap( M1 * M2   );  	 
' )

cat('Compliling calc_M1_M2_M3_Hadamard\n')
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

cat('Compliling calc_y_a\n')
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


cat('Compliling calc_tapply_vect_sum\n')
calc_tapply_vect_sum <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(v_i1 = "array", v_i2 = "array"), body='	
	// Calculate the equivalent of tapply(v_i1, v_i2, sum)
 	using Eigen::Map;		 		// to map input variable to an existing array of data 
	using Eigen::ArrayXd; 				// to use ArrayXd 
	using Eigen::ArrayXi; 				// to use ArrayXi (integers)  
  // eta.sn <- as.vector(Wtime2 %*% phi.new) + alpha.new * Ztime2.b # M*GQ matrix #
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


cat('Compliling calc_VY\n')
calc_VY <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(v_i1 = "matrix", a_i = "matrix", b_i = "numeric"), body='	
	// Calculate VY as $v_i * a_i * v_i^T + I * b_i$ 
 	using Eigen::Map;		 		// to map input variable to an existing array of data 
	using Eigen::ArrayXd; 				// to use ArrayXd  
	using Eigen::MatrixXd; 				// to use MatrixXd  

	const 	Map<MatrixXd> v(Rcpp::as<Map<MatrixXd> > (v_i1));	//Map array v_i1 to ArrayXd v1  
	const 	Map<MatrixXd> a(Rcpp::as<Map<MatrixXd> > (a_i));	//Map double a_i to double a 
	const double b (Rcpp::as<double> (b_i));			//Map double b_i to double b 

	MatrixXd M1 = v*a*v.transpose();
	M1.diagonal().array() += b;

	return Rcpp::wrap( M1 );  
' )


cat('Compliling calc_VB\n')
calc_VB <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 = "matrix", M_i2 ="matrix", M_i3 = "matrix"), body='
	// Calculate VB  as   a_i  -  a_i * y_i^T * M_i^(-1) * y_i * a_i
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition 

	const Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));	// Map matrix M_i to matrixXd M  
	const Map<MatrixXd> M2(Rcpp::as<Map<MatrixXd> > (M_i2));	// Map matrix M_i to matrixXd M  
	const Map<MatrixXd> M3(Rcpp::as<Map<MatrixXd> > (M_i3));	// Map matrix M_i to matrixXd M  
 	LLT<MatrixXd> llt_M3(M3);					// Calculate the LLT decomposition 
	
	MatrixXd b = M2.transpose() * llt_M3.solve(M2) ; //return Rcpp::wrap( b  );
	return Rcpp::wrap( M1 - M1 * b * M1   );  
	 
' )



cat('Compliling calc_muB\n')
calc_muB <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(  M_i3 = "matrix", y_i1 = "vector", y_i2 = "vector", M_i0 ="matrix", M_i1 ="matrix", M_i2 ="matrix"), body='
	//  as.vector(BSigma.old %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - X.st[[i]] %*% beta.old)))
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition   
 
	const Map<MatrixXd> BSold(Rcpp::as<Map<MatrixXd> > (M_i0));	// Map matrix M_i0 to matrixXd M
	const Map<MatrixXd> VY(Rcpp::as<Map<MatrixXd> > (M_i1));	// Map matrix M_i1 to matrixXd M
	const Map<MatrixXd> Xst(Rcpp::as<Map<MatrixXd> > (M_i2));	// Map matrix M_i2 to matrixXd M
	const Map<MatrixXd> Zst(Rcpp::as<Map<MatrixXd> > (M_i3));	// Map matrix M_i3 to vectorXd M  
	const Map<VectorXd> Yst(Rcpp::as<Map<VectorXd> > (y_i1));	// Map vector y_i1 to vectorXd y  
	const Map<VectorXd> betaold(Rcpp::as<Map<VectorXd> > (y_i2));	// Map vector y_i2 to vectorXd y  

	VectorXd yf = Yst - Xst * betaold;			
 	LLT<MatrixXd> llt_VY(VY); 					// Calculate the LLT decomposition 
	    	 
 	return Rcpp::wrap( BSold * Zst.transpose() * llt_VY.solve(yf) );  	 
' ) 
 


cat('Compliling calc_bi_st\n')
calc_bi_st <- cxxfunction( plugin = "RcppEigen",  signature(y_i0 = "vector", y_i1 = "matrix", M_i ="matrix"), body='
	// Calculate $chol(M_i^{-1})$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::LLT; 				// to do the LLT decomposition  
	using Eigen::Lower;
  
	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i));		// Map matrix M_i to matrixXd M
	const Map<MatrixXd> y1(Rcpp::as<Map<MatrixXd> > (y_i1));	// Map vector y_i to vectorXd y  
	const Map<VectorXd> y0(Rcpp::as<Map<VectorXd> > (y_i0));	// Map vector y_i to vectorXd y  
	
	LLT<MatrixXd> llt_M(M.inverse()); 
	MatrixXd Res = sqrt(2.) *llt_M.matrixLLT().triangularView<Lower>().transpose().solve(y1.transpose());

 	for (unsigned int i = 0; i != Res.cols(); ++i) {
		Res.col(i) += y0;		
	}

	return Rcpp::wrap( Res  );
' )

cat('Compliling calc_MVND\n')
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

cat('Compliling calc_rowsum\n')
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

cat('Compliling calc_mult_rowsum\n')
calc_mult_rowsum <- cxxfunction(settings=settingsE, plugin="RcppEigen", signature(y_i="vector", y_i2 = "vector", M_i1="Matrix", M_i2="Array"), body='
	//CondExp2 * calc_rowsum( y_i = Index , Xtime22[, i] * temp0d))
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXi; 				// to use VectorXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::ArrayXd;
	//#include <iostream>

	const Map<VectorXi> y1(Rcpp::as<Map<VectorXi> > (y_i));		// Map vector y_i1 to vectorXd y1  
	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i1));		// Map matrix M_i1 to vectorXd M  
	const Map<ArrayXd> A(Rcpp::as<Map<ArrayXd> > (M_i2));		// Map matrix M_i2 to vectorXd M  
	const Map<VectorXd> y2(Rcpp::as<Map<VectorXd> > (y_i2));	// Map vector y_i1 to vectorXd y1  
	
	const unsigned int l = y1.size();
	const unsigned int m = M.cols(); 	

	MatrixXd Res = MatrixXd::Zero(y1.maxCoeff() ,m);	

	for (unsigned int i = 0; i < m; ++i){
	unsigned int k = 0;
		for(unsigned int j = 0; j != l; ++j){
 			Res(k,i) = Res(k,i) + M(j,i) * y2(j); 
			if (j  < M.rows() ) {
				if( y1[j] != y1[j+1] ){ 
					k++;
				} 
			}
		}
	}

	Res.array() *= A;

	return Rcpp::wrap( Res );  	 
' )


 
cat('Compliling calc_mult_rowsum2\n')
calc_mult_rowsum2 <- cxxfunction(settings=settingsE, plugin="RcppEigen", signature(y_i="vector", M_i3 = "matrix", M_i1="Matrix", M_i2="Array"), body='
	//CondExp2 * calc_rowsum( y_i = Index , Xtime22 * temp0d))
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXi; 				// to use VectorXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::ArrayXd;
	//#include <iostream>

	const Map<VectorXi> y1(Rcpp::as<Map<VectorXi> > (y_i));		// Map vector y_i1 to vectorXd y1  
	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i1));		// Map matrix M_i1 to vectorXd M  
	const Map<ArrayXd> A(Rcpp::as<Map<ArrayXd> > (M_i2));		// Map matrix M_i2 to vectorXd M  
	const Map<MatrixXd> M3(Rcpp::as<Map<MatrixXd> > (M_i3));	// Map vector y_i1 to vectorXd y1  
	
	const unsigned int l = y1.size();
	const unsigned int m = M.cols(); 	

	MatrixXd Res = MatrixXd::Zero(y1.maxCoeff() ,m);	

	for (unsigned int i = 0; i < m; ++i){
	unsigned int k = 0;
		for(unsigned int j = 0; j != l; ++j){
 			Res(k,i) = Res(k,i) + M(j,i) * M3(j,i); 
			if (j  < M.rows() ) {
				if( y1[j] != y1[j+1] ){ 
					k++;
				} 
			}
		}
	}

	Res.array() *= A;

	return Rcpp::wrap( Res );  	 
' )



cat('Compliling calc_mult0_rowsum\n')
calc_mult0_rowsum <- cxxfunction(settings=settingsE, plugin="RcppEigen", signature(y_i="vector", y_i2 = "vector", M_i1="Matrix" ), body='
	//  calc_rowsum( y_i = Index , Xtime22[, i] * temp0d))
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXi; 				// to use VectorXd
	using Eigen::VectorXd; 				// to use VectorXd
	using Eigen::ArrayXd;
	//#include <iostream>

	const Map<VectorXi> y1(Rcpp::as<Map<VectorXi> > (y_i));		// Map vector y_i1 to vectorXd y1  
	const Map<MatrixXd> M(Rcpp::as<Map<MatrixXd> > (M_i1));		// Map matrix M_i1 to vectorXd M   
	const Map<VectorXd> y2(Rcpp::as<Map<VectorXd> > (y_i2));	// Map vector y_i1 to vectorXd y1  
	
	const unsigned int l = y1.size();
	const unsigned int m = M.cols(); 	

	MatrixXd Res = MatrixXd::Zero(y1.maxCoeff() ,m);	

	for (unsigned int i = 0; i < m; ++i){
	unsigned int k = 0;
		for(unsigned int j = 0; j != l; ++j){
 			Res(k,i) = Res(k,i) + M(j,i) * y2(j); 
			if (j  < M.rows() ) {
				if( y1[j] != y1[j+1] ){ 
					k++;
				} 
			}
		}
	}
 

	return Rcpp::wrap( Res );  	 
' )


 src <- '
using Eigen::MatrixXd; 				// to use MatrixXd
using Eigen::VectorXd; 				// to use MatrixXd
using Eigen::Map;

Rcpp::List const input(data); 

const unsigned int l = input.size();

Rcpp::NumericMatrix xx = input[1];

const unsigned int nr = xx.nrow();
const unsigned int nc = xx.ncol();

MatrixXd U(nr*nr*l, nc);
MatrixXd u1 =  MatrixXd::Zero( nr, nc );	
MatrixXd u2 =  MatrixXd::Zero( nr, nr );	

for (unsigned int i = 0; i != l; ++i){  
  u1 =  input[i];
  for (unsigned int j = 0; j != nc; ++j){
    u2 = u1.col(j) *  u1.col(j).adjoint();
    U.block( i*nr*nr, j, nr*nr, 1) = VectorXd::Map(u2.data(), u2.rows()*u2.cols()); 
   //U.block( i*nr, j, nr*nr, 1) = VectorXd::Map(u2.data(), u2.rows()*u2.cols()); 
   // std::cout << u2 << std::endl;
 }  
} 
	return Rcpp::wrap( U );  	
'

cat('Compliling fast_rbind_lapply\n')
fast_rbind_lapply <- cxxfunction(signature(data = "list"),  src, plugin = "RcppEigen")


cat('Compliling calc_exp2M\n')
calc_expM2 <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i ="array"), body='	
	// Calculate $exp(M_i)$ - this function returns an array
 	using Eigen::Map;		 		// to map input variable to an existing array of data 
	using Eigen::ArrayXd; 				// to use ArrayXd 

	Map<ArrayXd> A(Rcpp::as<Map<ArrayXd> > (M_i));	//Map vector M_i to arrayXd A  
 	A = (A.exp());
	//return Rcpp::wrap(  );   
' )



 src <- '
using Eigen::MatrixXd; 				// to use MatrixXd
using Eigen::VectorXd; 				// to use MatrixXd
using Eigen::Map;

Rcpp::List const input1(data1); 
Rcpp::List const input2(data2); 
Rcpp::NumericVector const Ind(indeces);

const unsigned int l = Ind.size();
unsigned int s = 0; 

for (unsigned int i = 0; i != l; ++i){  
  Rcpp::NumericMatrix xx = input1[ Ind(i) ];
  s += xx.nrow();
}

Rcpp::NumericMatrix xx = input2[1];

const unsigned int nc = xx.ncol();
MatrixXd U( s , nc );
 
unsigned int j = 0;
for (unsigned int i = 0; i != l; ++i){  
	//const Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (input1[ Ind(i) ]));
	//const Map<MatrixXd> M2(Rcpp::as<Map<MatrixXd> > (input2[ Ind(i) ]));
	
	const MatrixXd M1 = (input1[ Ind(i) ]);
	const MatrixXd M2 = (input2[ Ind(i) ]);
	
	// std::cout << "M1.rows()" << M1.rows() << std::endl;
	// std::cout << "j = " <<j << std::endl;
	 U.block( j, 0, M1.rows() , nc ) =   M1 * M2;
 	j += M1.rows();
}
	// std::cout << "j = " <<j << std::endl;
	return Rcpp::wrap( U );  	
'

cat('Compliling fast_lapply_length \n')
fast_lapply_length <- cxxfunction(signature(data1 = "list", data2 = "list", indeces = "vector"),  src, plugin = "RcppEigen")


cat('Finished compiling\n')




