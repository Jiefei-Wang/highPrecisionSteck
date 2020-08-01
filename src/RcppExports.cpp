// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rational_steck
NumericVector rational_steck(NumericVector l, NumericVector h);
RcppExport SEXP _highPrecisionSteck_rational_steck(SEXP lSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(rational_steck(l, h));
    return rcpp_result_gen;
END_RCPP
}
// rational_row
NumericVector rational_row(NumericVector l, NumericVector h);
RcppExport SEXP _highPrecisionSteck_rational_row(SEXP lSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(rational_row(l, h));
    return rcpp_result_gen;
END_RCPP
}
// high_steck
NumericVector high_steck(NumericVector l, NumericVector h, int prec, NumericVector upperBound);
RcppExport SEXP _highPrecisionSteck_high_steck(SEXP lSEXP, SEXP hSEXP, SEXP precSEXP, SEXP upperBoundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type prec(precSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upperBound(upperBoundSEXP);
    rcpp_result_gen = Rcpp::wrap(high_steck(l, h, prec, upperBound));
    return rcpp_result_gen;
END_RCPP
}
// high_row
NumericVector high_row(NumericVector l, NumericVector h, int prec, NumericVector upperBound);
RcppExport SEXP _highPrecisionSteck_high_row(SEXP lSEXP, SEXP hSEXP, SEXP precSEXP, SEXP upperBoundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type prec(precSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upperBound(upperBoundSEXP);
    rcpp_result_gen = Rcpp::wrap(high_row(l, h, prec, upperBound));
    return rcpp_result_gen;
END_RCPP
}
// set_max_binomial
void set_max_binomial(int n);
RcppExport SEXP _highPrecisionSteck_set_max_binomial(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    set_max_binomial(n);
    return R_NilValue;
END_RCPP
}
// get_log_prod
std::vector<double> get_log_prod();
RcppExport SEXP _highPrecisionSteck_get_log_prod() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(get_log_prod());
    return rcpp_result_gen;
END_RCPP
}
// C_log_prod
double C_log_prod(int n);
RcppExport SEXP _highPrecisionSteck_C_log_prod(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(C_log_prod(n));
    return rcpp_result_gen;
END_RCPP
}
// C_log_beta
double C_log_beta(int a, int b);
RcppExport SEXP _highPrecisionSteck_C_log_beta(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(C_log_beta(a, b));
    return rcpp_result_gen;
END_RCPP
}
// C_lchoose
double C_lchoose(int n, int k);
RcppExport SEXP _highPrecisionSteck_C_lchoose(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(C_lchoose(n, k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_highPrecisionSteck_rational_steck", (DL_FUNC) &_highPrecisionSteck_rational_steck, 2},
    {"_highPrecisionSteck_rational_row", (DL_FUNC) &_highPrecisionSteck_rational_row, 2},
    {"_highPrecisionSteck_high_steck", (DL_FUNC) &_highPrecisionSteck_high_steck, 4},
    {"_highPrecisionSteck_high_row", (DL_FUNC) &_highPrecisionSteck_high_row, 4},
    {"_highPrecisionSteck_set_max_binomial", (DL_FUNC) &_highPrecisionSteck_set_max_binomial, 1},
    {"_highPrecisionSteck_get_log_prod", (DL_FUNC) &_highPrecisionSteck_get_log_prod, 0},
    {"_highPrecisionSteck_C_log_prod", (DL_FUNC) &_highPrecisionSteck_C_log_prod, 1},
    {"_highPrecisionSteck_C_log_beta", (DL_FUNC) &_highPrecisionSteck_C_log_beta, 2},
    {"_highPrecisionSteck_C_lchoose", (DL_FUNC) &_highPrecisionSteck_C_lchoose, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_highPrecisionSteck(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
