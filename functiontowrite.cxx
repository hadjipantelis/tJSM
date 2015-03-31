
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
}

if(1==3){
# // tempB <- do.call(rbind, lapply(1:n, function(i) apply((bi.st[[i]]), 2, function(x) tcrossprod(x) )))
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

cpp_testy <- cxxfunction(signature(data = "list"),  src, plugin = "RcppEigen")


cat('Compliling calc_M_M\n')
calc_M_M <- cxxfunction(settings=settingsE, plugin = "RcppEigen", signature(M_i1 = "matrix", M_i2 ="matrix"), body='	
	// Calculate $M_i y_i$
 	using Eigen::Map;		 		// to map input variable to an existing array of data
	using Eigen::MatrixXd; 				// to use MatrixXd 

	const Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));	//Map vector M_i to MatrixXd M 
	const Map<MatrixXd> M2(Rcpp::as<Map<MatrixXd> > (M_i2));	//Map vector y_i to vectorXd y 

	return Rcpp::wrap(  M1 * M2  );  	 
' )



cat('Compliling calc_etan\n')
calc_etan <- cxxfunction(settings=settingsE, plugin = "RcppEigen",  signature(M_i1 = "matrix", M_i2 ="array", a_i1 = "numeric", y_i = "vector"), body='	
	// Calculate $exp(M_i)$ - this function returns an array
 	using Eigen::Map;		 		// to map input variable to an existing array of data 
	using Eigen::ArrayXXd; 				// to use ArrayXd 
	using Eigen::MatrixXd; 				// to use ArrayXd 

  // eta.sn <- as.vector(Wtime2 %*% phi.new) + alpha.new * Ztime2.b # M*GQ matrix #
	const double a1 (Rcpp::as<double> (a_i1));		//Map double a_i to double a 
	const Map<MatrixXd> M1(Rcpp::as<Map<MatrixXd> > (M_i1));//Map vector M_i to MatrixXd M 
		Map<ArrayXXd> M2(Rcpp::as<Map<ArrayXXd> > (M_i2));//Map vector y_i to vectorXd y 
	const Map<MatrixXd> y(Rcpp::as<Map<MatrixXd> > (y_i));	//Map vector y_i to vectorXd y 
  
	MatrixXd U = M1 * y;
	M2*= a1;

	for (unsigned int k=0; k!= M2.rows(); ++k){

	M2.row(k)  += U(k);
}  
' )

}
 

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

