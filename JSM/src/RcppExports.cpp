// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// calc_bi_st
Eigen::MatrixXd calc_bi_st(const Eigen::Map<Eigen::VectorXd>& v0, const Eigen::Map<Eigen::MatrixXd>& m, const Eigen::Map<Eigen::MatrixXd>& M);
RcppExport SEXP JSM_calc_bi_st(SEXP v0SEXP, SEXP mSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M(MSEXP);
    __result = Rcpp::wrap(calc_bi_st(v0, m, M));
    return __result;
END_RCPP
}
// calc_expM2
void calc_expM2(Eigen::Map<Eigen::ArrayXd>& A);
RcppExport SEXP JSM_calc_expM2(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd>& >::type A(ASEXP);
    calc_expM2(A);
    return R_NilValue;
END_RCPP
}
// calc_M1_a_M2_Hadamard
void calc_M1_a_M2_Hadamard(Eigen::Map<Eigen::MatrixXd>& M1, const Eigen::Map<Eigen::MatrixXd>& M2, const double a, const Eigen::Map<Eigen::VectorXi>& v);
RcppExport SEXP JSM_calc_M1_a_M2_Hadamard(SEXP M1SEXP, SEXP M2SEXP, SEXP aSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M2(M2SEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXi>& >::type v(vSEXP);
    calc_M1_a_M2_Hadamard(M1, M2, a, v);
    return R_NilValue;
END_RCPP
}
// calc_M1_M2_Hadamard_a
Eigen::MatrixXd calc_M1_M2_Hadamard_a(Eigen::Map<Eigen::ArrayXXd>& A1, const Eigen::Map<Eigen::ArrayXXd>& A2, const Eigen::Map<Eigen::VectorXd>& v3, const int a);
RcppExport SEXP JSM_calc_M1_M2_Hadamard_a(SEXP A1SEXP, SEXP A2SEXP, SEXP v3SEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXXd>& >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXXd>& >::type A2(A2SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type v3(v3SEXP);
    Rcpp::traits::input_parameter< const int >::type a(aSEXP);
    __result = Rcpp::wrap(calc_M1_M2_Hadamard_a(A1, A2, v3, a));
    return __result;
END_RCPP
}
// calc_M1_M2_Hadamard
void calc_M1_M2_Hadamard(Eigen::Map<Eigen::ArrayXd>& M1, const Eigen::Map<Eigen::ArrayXd>& M2);
RcppExport SEXP JSM_calc_M1_M2_Hadamard(SEXP M1SEXP, SEXP M2SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd>& >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd>& >::type M2(M2SEXP);
    calc_M1_M2_Hadamard(M1, M2);
    return R_NilValue;
END_RCPP
}
// calc_M1_M2_M3_Hadamard
void calc_M1_M2_M3_Hadamard(Eigen::Map<Eigen::MatrixXd>& M1, const Eigen::Map<Eigen::MatrixXd>& M2, const Eigen::Map<Eigen::MatrixXd>& M3, const Eigen::Map<Eigen::VectorXi>& v);
RcppExport SEXP JSM_calc_M1_M2_M3_Hadamard(SEXP M1SEXP, SEXP M2SEXP, SEXP M3SEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M2(M2SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M3(M3SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXi>& >::type v(vSEXP);
    calc_M1_M2_M3_Hadamard(M1, M2, M3, v);
    return R_NilValue;
END_RCPP
}
// calc_M1timesM2v
Eigen::MatrixXd calc_M1timesM2v(const Eigen::Map<Eigen::MatrixXd>& M1, const Eigen::Map<Eigen::MatrixXd>& M2, const Eigen::Map<Eigen::ArrayXd>& v);
RcppExport SEXP JSM_calc_M1timesM2v(SEXP M1SEXP, SEXP M2SEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M2(M2SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd>& >::type v(vSEXP);
    __result = Rcpp::wrap(calc_M1timesM2v(M1, M2, v));
    return __result;
END_RCPP
}
// calc_muB
Eigen::MatrixXd calc_muB(const Eigen::Map<Eigen::MatrixXd>& BSold, const Eigen::Map<Eigen::MatrixXd>& VY, const Eigen::Map<Eigen::MatrixXd>& Xst, const Eigen::Map<Eigen::MatrixXd>& Zst, const Eigen::Map<Eigen::VectorXd>& Yst, const Eigen::Map<Eigen::VectorXd>& betaold);
RcppExport SEXP JSM_calc_muB(SEXP BSoldSEXP, SEXP VYSEXP, SEXP XstSEXP, SEXP ZstSEXP, SEXP YstSEXP, SEXP betaoldSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type BSold(BSoldSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type VY(VYSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Xst(XstSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Zst(ZstSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type Yst(YstSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type betaold(betaoldSEXP);
    __result = Rcpp::wrap(calc_muB(BSold, VY, Xst, Zst, Yst, betaold));
    return __result;
END_RCPP
}
// calc_muBMult
Eigen::VectorXd calc_muBMult(const Eigen::Map<Eigen::MatrixXd>& BSold, const Eigen::Map<Eigen::MatrixXd>& VY, const Eigen::Map<Eigen::VectorXd>& BTg, const Eigen::Map<Eigen::VectorXd>& Yst);
RcppExport SEXP JSM_calc_muBMult(SEXP BSoldSEXP, SEXP VYSEXP, SEXP BTgSEXP, SEXP YstSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type BSold(BSoldSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type VY(VYSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type BTg(BTgSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type Yst(YstSEXP);
    __result = Rcpp::wrap(calc_muBMult(BSold, VY, BTg, Yst));
    return __result;
END_RCPP
}
// calc_mult_rowsum1
Eigen::MatrixXd calc_mult_rowsum1(const Eigen::Map<Eigen::VectorXi>& v, const Eigen::Map<Eigen::VectorXd>& u, const Eigen::Map<Eigen::MatrixXd>& M, const Eigen::Map<Eigen::ArrayXd>& A);
RcppExport SEXP JSM_calc_mult_rowsum1(SEXP vSEXP, SEXP uSEXP, SEXP MSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXi>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd>& >::type A(ASEXP);
    __result = Rcpp::wrap(calc_mult_rowsum1(v, u, M, A));
    return __result;
END_RCPP
}
// calc_mult_rowsum2
Eigen::MatrixXd calc_mult_rowsum2(const Eigen::Map<Eigen::VectorXi>& v, const Eigen::Map<Eigen::MatrixXd>& L, const Eigen::Map<Eigen::MatrixXd>& M, const Eigen::Map<Eigen::ArrayXd>& A);
RcppExport SEXP JSM_calc_mult_rowsum2(SEXP vSEXP, SEXP LSEXP, SEXP MSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXi>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd>& >::type A(ASEXP);
    __result = Rcpp::wrap(calc_mult_rowsum2(v, L, M, A));
    return __result;
END_RCPP
}
// calc_mult_rowsum3
Rcpp::List calc_mult_rowsum3(const Eigen::Map<Eigen::ArrayXi>& v, const Eigen::Map<Eigen::ArrayXXd>& B, const Eigen::Map<Eigen::ArrayXXd>& M, const Eigen::Map<Eigen::ArrayXd>& A, const double ncb2);
RcppExport SEXP JSM_calc_mult_rowsum3(SEXP vSEXP, SEXP BSEXP, SEXP MSEXP, SEXP ASEXP, SEXP ncb2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXi>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXXd>& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXXd>& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd>& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const double >::type ncb2(ncb2SEXP);
    __result = Rcpp::wrap(calc_mult_rowsum3(v, B, M, A, ncb2));
    return __result;
END_RCPP
}
// calc_M_v
Eigen::VectorXd calc_M_v(const Eigen::Map<Eigen::VectorXd>& v, const Eigen::Map<Eigen::MatrixXd>& M);
RcppExport SEXP JSM_calc_M_v(SEXP vSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M(MSEXP);
    __result = Rcpp::wrap(calc_M_v(v, M));
    return __result;
END_RCPP
}
// calc_MVND
double calc_MVND(const Eigen::Map<Eigen::VectorXd>& x, const Eigen::Map<Eigen::VectorXd>& mu, const Eigen::Map<Eigen::MatrixXd>& K);
RcppExport SEXP JSM_calc_MVND(SEXP xSEXP, SEXP muSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type K(KSEXP);
    __result = Rcpp::wrap(calc_MVND(x, mu, K));
    return __result;
END_RCPP
}
// calc_rowsum
Eigen::MatrixXd calc_rowsum(const Eigen::Map<Eigen::VectorXi>& v, const Eigen::Map<Eigen::MatrixXd>& M);
RcppExport SEXP JSM_calc_rowsum(SEXP vSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXi>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M(MSEXP);
    __result = Rcpp::wrap(calc_rowsum(v, M));
    return __result;
END_RCPP
}
// calc_rowsum_mult
Eigen::MatrixXd calc_rowsum_mult(const Eigen::Map<Eigen::VectorXi>& v, const Eigen::Map<Eigen::VectorXd>& u, const Eigen::Map<Eigen::MatrixXd>& M);
RcppExport SEXP JSM_calc_rowsum_mult(SEXP vSEXP, SEXP uSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXi>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M(MSEXP);
    __result = Rcpp::wrap(calc_rowsum_mult(v, u, M));
    return __result;
END_RCPP
}
// calc_tapply_vect_sum
Eigen::ArrayXd calc_tapply_vect_sum(const Eigen::Map<Eigen::ArrayXd>& v1, const Eigen::Map<Eigen::ArrayXi>& v2);
RcppExport SEXP JSM_calc_tapply_vect_sum(SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXd>& >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::ArrayXi>& >::type v2(v2SEXP);
    __result = Rcpp::wrap(calc_tapply_vect_sum(v1, v2));
    return __result;
END_RCPP
}
// calc_v_a
void calc_v_a(Eigen::Map<Eigen::ArrayXd>& v, const double& a);
RcppExport SEXP JSM_calc_v_a(SEXP vSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd>& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    calc_v_a(v, a);
    return R_NilValue;
END_RCPP
}
// calc_VB
Eigen::MatrixXd calc_VB(const Eigen::Map<Eigen::MatrixXd>& M1, const Eigen::Map<Eigen::MatrixXd>& M2, const Eigen::Map<Eigen::MatrixXd>& M3);
RcppExport SEXP JSM_calc_VB(SEXP M1SEXP, SEXP M2SEXP, SEXP M3SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M2(M2SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M3(M3SEXP);
    __result = Rcpp::wrap(calc_VB(M1, M2, M3));
    return __result;
END_RCPP
}
// calc_VY
Eigen::MatrixXd calc_VY(const Eigen::Map<Eigen::MatrixXd>& M, const Eigen::Map<Eigen::MatrixXd>& A, const double b);
RcppExport SEXP JSM_calc_VY(SEXP MSEXP, SEXP ASEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    __result = Rcpp::wrap(calc_VY(M, A, b));
    return __result;
END_RCPP
}
// fast_lapply_length
Eigen::MatrixXd fast_lapply_length(Rcpp::List const input1, Rcpp::List const input2, Rcpp::NumericVector const Ind);
RcppExport SEXP JSM_fast_lapply_length(SEXP input1SEXP, SEXP input2SEXP, SEXP IndSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List const >::type input1(input1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List const >::type input2(input2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const >::type Ind(IndSEXP);
    __result = Rcpp::wrap(fast_lapply_length(input1, input2, Ind));
    return __result;
END_RCPP
}
// fast_rbind_lapply_outerprod
Eigen::MatrixXd fast_rbind_lapply_outerprod(Rcpp::List const input);
RcppExport SEXP JSM_fast_rbind_lapply_outerprod(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List const >::type input(inputSEXP);
    __result = Rcpp::wrap(fast_rbind_lapply_outerprod(input));
    return __result;
END_RCPP
}
