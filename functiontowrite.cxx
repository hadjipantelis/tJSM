
if(11==1){
print('Compliling calc_bi_st\n')
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
if(1==11){
print('Compliling calc_bi_st\n')
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

}

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
