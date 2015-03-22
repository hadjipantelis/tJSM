
if(11==1){
  
 
cat('Compliling calc_tcrossprod\n')
calc_tcrossprod <- cxxfunction( plugin = "RcppEigen",  signature( M_i0 ="matrix"), body='
	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::Lower; 				// to use  
 
	const Map<MatrixXd> A(Rcpp::as<Map<MatrixXd> > (M_i0));	// Map matrix M_i0 to matrixXd M
	const unsigned int m(A.rows());	 

	MatrixXd Res(MatrixXd(m,m).setZero().selfadjointView<Lower>().rankUpdate(A));
 	return Rcpp::wrap(Res);  	 
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
}



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




