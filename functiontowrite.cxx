print('Compliling calc_rowsum')
calc_rowsum <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(y_i = "vector", M_i = "Matrix"), body='
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd
	using Eigen::VectorXd; 				// to use VectorXd
	//#include <iostream>

	const Map<VectorXd> y1(Rcpp::as<Map<VectorXd> > (y_i));	// Map vector y_i1 to vectorXd y1  
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


rowsum_f <- function( A,b){

	l = length(b);
	M = ncol(A);
	Res = matrix(rep(0, M * max(b)), ncol= M)

	for (i in 1:M){
	k = 1;
 		for (j in 1:(l)) {
        		Res[k,i] = Res[k,i] + A[j,i];
			if ( j +1  <= l ) {
				if ( b[j] != b[j+1]) { 
					k = k+1
				} 
			}
		}
	} 
	return (Res)
}
