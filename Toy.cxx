
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

code='
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::export]]
double long_computation_omp(int nb, int threads=1) {
#ifdef _OPENMP
    if ( threads > 0 )
        omp_set_num_threads( threads );
    REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif
 
    double sum = 0;
#pragma omp parallel for schedule(dynamic)   
    for (int i = 0; i < nb; ++i) {
        double thread_sum = 0;
  	for (int j = 0; j < nb; ++j) {
	    thread_sum += R::dlnorm(i+j, 0.0, 1.0, 0);
	}
        sum += thread_sum;
    }
    return sum + nb;
}
'

sourceCpp(code=code)
s <- long_computation_omp2(10000, 4)

settings <- getPlugin("RcppEigen")
settings$env$PKG_CXXFLAGS <- paste('-fopenmp', settings$env$PKG_CXXFLAGS)
settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)


 
print('Compliling fun_par')
src <-
' 
#ifdef _OPENMP
#include <omp.h>
#endif
Rcpp::NumericVector xa(a);
Rcpp::NumericVector xb(b);
Rcpp::NumericVector xc(c);
int n_xa = xa.size(), n_xb = xb.size(), n_xc = xc.size();  
#pragma omp parallel for schedule(dynamic)   
for (int i = 0; i < n_xa; i++)
xa[i] *= xb[i] * xc[i];
return xa; 
'

fun_par <- cxxfunction(signature(a = "numeric", b = "numeric", c= "numeric"), src, plugin = "Rcpp")


print('Compliling calc_M1_M2_M3')
calc_M1_M2_M3 <- cxxfunction(settings=settings, plugin = "RcppEigen",  signature(M_i1 ="matrix", M_i2 ="matrix", M_i3 ="matrix"), body='	
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


