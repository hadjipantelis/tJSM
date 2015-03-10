if(1==12){
print('Compliling calc_M1_M2_M3_b')
calc_M1_M2_M3_b <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 ="matrix", M_i2 ="matrix", M_i3 ="matrix", v_i = "vector"), body='	
	// Calculate the product $M1 M2 M3$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd 
	using Eigen::VectorXi; 				// to use VectorXd 

	Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));	//Map matrix M_1 to MatrixXd M1
	const Map<MatrixXd> M2(Rcpp::as<Map<MatrixXd> > (M_i2));//Map matrix M_2 to MatrixXd M2
	const Map<MatrixXd> M3(Rcpp::as<Map<MatrixXd> > (M_i3));//Map matrix M_3 to MatrixXd M3
	const Map<VectorXi> v(Rcpp::as<Map<VectorXi> > (v_i));	//Map vector v_i to MatrixXd v
	int N = v.size();

	for(unsigned int i=0; i != N; ++i){
		 M1.row(i).array()  *= M2.row(v(i)).array() *  M3.row(v(i)).array();
	}
 
	//return Rcpp::wrap( M1 );  	 
' )

if(1==1){
print('Compliling calc_M1_M2_M3_b')
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

}

print('Compliling calc_M1_M2_M3_c')
calc_M1_M2_M3_c <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 ="array", M_i2 ="array", M_i3 ="array", v_i = "vector", M = "numeric"), body='	
	// Calculate the product $M1 M2 M3$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::ArrayXd; 				// to use MatrixXd 
	using Eigen::VectorXi; 				// to use VectorXd 
	#include <iostream> 
	      Map<ArrayXd> M1(Rcpp::as<Map<ArrayXd> > (M_i1));	//Map matrix M_1 to MatrixXd M1
	const Map<ArrayXd> M2(Rcpp::as<Map<ArrayXd> > (M_i2));//Map matrix M_2 to MatrixXd M2
	const Map<ArrayXd> M3(Rcpp::as<Map<ArrayXd> > (M_i3));//Map matrix M_3 to MatrixXd M3
	const Map<VectorXi> v(Rcpp::as<Map<VectorXi> > (v_i));	//Map integer to integer
	const int m (Rcpp::as<int> (M));	//Map vector v_i to MatrixXd v
	int N = v.size();

	for(unsigned int i=0; i < N; i=m+i){
		 M1.segment(i,m) *= M2.segment(i,m) * M3.segment(i,m)  ;
	}
 
	//return Rcpp::wrap( M1  );  	 
' )


}
calc_fun <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(v_i1 ="vector", v_i2 ="vector", v_i3 ="vector"), body='	
	// Calculate the product $M1 M2 M3$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
 	using Eigen::VectorXi; 				// to use VectorXI 
	using Eigen::VectorXd; 				// to use VectorXd 
	const Map<VectorXd> v1(Rcpp::as<Map<VectorXd> > (v_i1));	//Map integer to integer
	const Map<VectorXi> v2(Rcpp::as<Map<VectorXi> > (v_i2));	//Map integer to integer
	const Map<VectorXi> v3(Rcpp::as<Map<VectorXi> > (v_i3));	//Map integer to integer
	double temp =0; 
	VectorXd Z = VectorXd::Zero( v3.size());
	for(unsigned int i=0; i <  v3.size(); i=1+i){ 
		for (unsigned int j=0; j <  v2.size(); j=1+i){
			if ( v3(i) == v2(j) ){
				temp += v1(j);			
			}				
		}
		Z(i) = temp; temp=0;
	}
	
	return Rcpp::wrap( Z );
' )
 

#calc_M1_M2_M3_c( exp.es2, CondExp , Integral,as.integer(Index-1),12);print(exp.es2[1:3]);exp.es2 <- exp.es; exp.es2[1]= exp.es2[1] +0






