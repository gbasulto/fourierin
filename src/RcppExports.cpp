// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ecf_mod
NumericVector ecf_mod(NumericVector t, NumericVector smp);
RcppExport SEXP fourierin_ecf_mod(SEXP tSEXP, SEXP smpSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type smp(smpSEXP);
    __result = Rcpp::wrap(ecf_mod(t, smp));
    return __result;
END_RCPP
}
