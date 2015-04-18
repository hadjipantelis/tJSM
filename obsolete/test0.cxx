
print('Compliling calc_M1_M2_M3')
calc_M1_M2_M3 <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 ="matrix", M_i2 ="matrix", M_i3 ="matrix"), body='	
	// Calculate the product $M1 M2 M3$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd 
	#include <iostream> 
	const Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));	//Map vector M_1 to MatrixXd M1
	const Map<MatrixXd> M2(Rcpp::as<Map<MatrixXd> > (M_i2));	//Map vector M_2 to MatrixXd M2
	const Map<MatrixXd> M3(Rcpp::as<Map<MatrixXd> > (M_i3));	//Map vector M_3 to MatrixXd M3 
	std::cout<< Eigen::nbThreads() << std::endl;
	return Rcpp::wrap( M1.cwiseProduct(M2).cwiseProduct(M3)  );  	 
' )