// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/fourierin.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// fourierin_1d_cpp
arma::cx_vec fourierin_1d_cpp(const arma::vec& f, double a, double b, double c, double d, double r, double s);
static SEXP fourierin_fourierin_1d_cpp_try(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP rSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const arma::vec& >::type f(fSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    __result = Rcpp::wrap(fourierin_1d_cpp(f, a, b, c, d, r, s));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP fourierin_fourierin_1d_cpp(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP rSEXP, SEXP sSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(fourierin_fourierin_1d_cpp_try(fSEXP, aSEXP, bSEXP, cSEXP, dSEXP, rSEXP, sSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// fourierin_1d_nonregular_cpp
arma::cx_vec fourierin_1d_nonregular_cpp(const arma::vec& f, double a, double b, const arma::vec& w, int resolution, double r, double s);
static SEXP fourierin_fourierin_1d_nonregular_cpp_try(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP wSEXP, SEXP resolutionSEXP, SEXP rSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const arma::vec& >::type f(fSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type resolution(resolutionSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    __result = Rcpp::wrap(fourierin_1d_nonregular_cpp(f, a, b, w, resolution, r, s));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP fourierin_fourierin_1d_nonregular_cpp(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP wSEXP, SEXP resolutionSEXP, SEXP rSEXP, SEXP sSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(fourierin_fourierin_1d_nonregular_cpp_try(fSEXP, aSEXP, bSEXP, wSEXP, resolutionSEXP, rSEXP, sSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// fourierin_2d_cpp
arma::cx_mat fourierin_2d_cpp(const arma::mat& f, const arma::vec& a, const arma::vec& b, const arma::vec& c, const arma::vec& d, double r, double s);
static SEXP fourierin_fourierin_2d_cpp_try(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP rSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const arma::mat& >::type f(fSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    __result = Rcpp::wrap(fourierin_2d_cpp(f, a, b, c, d, r, s));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP fourierin_fourierin_2d_cpp(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP rSEXP, SEXP sSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(fourierin_fourierin_2d_cpp_try(fSEXP, aSEXP, bSEXP, cSEXP, dSEXP, rSEXP, sSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// fourierin_cx_1d_cpp
arma::cx_vec fourierin_cx_1d_cpp(const arma::cx_vec& f, double a, double b, double c, double d, double r, double s);
static SEXP fourierin_fourierin_cx_1d_cpp_try(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP rSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type f(fSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    __result = Rcpp::wrap(fourierin_cx_1d_cpp(f, a, b, c, d, r, s));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP fourierin_fourierin_cx_1d_cpp(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP rSEXP, SEXP sSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(fourierin_fourierin_cx_1d_cpp_try(fSEXP, aSEXP, bSEXP, cSEXP, dSEXP, rSEXP, sSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// fourierin_cx_1d_nonregular_cpp
arma::cx_vec fourierin_cx_1d_nonregular_cpp(const arma::cx_vec& f, double a, double b, const arma::vec& w, int resolution, double r, double s);
static SEXP fourierin_fourierin_cx_1d_nonregular_cpp_try(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP wSEXP, SEXP resolutionSEXP, SEXP rSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type f(fSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type resolution(resolutionSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    __result = Rcpp::wrap(fourierin_cx_1d_nonregular_cpp(f, a, b, w, resolution, r, s));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP fourierin_fourierin_cx_1d_nonregular_cpp(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP wSEXP, SEXP resolutionSEXP, SEXP rSEXP, SEXP sSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(fourierin_fourierin_cx_1d_nonregular_cpp_try(fSEXP, aSEXP, bSEXP, wSEXP, resolutionSEXP, rSEXP, sSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// fourierin_cx_2d_cpp
arma::cx_mat fourierin_cx_2d_cpp(const arma::cx_mat& f, const arma::vec& a, const arma::vec& b, const arma::vec& c, const arma::vec& d, double r, double s);
static SEXP fourierin_fourierin_cx_2d_cpp_try(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP rSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type f(fSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    __result = Rcpp::wrap(fourierin_cx_2d_cpp(f, a, b, c, d, r, s));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP fourierin_fourierin_cx_2d_cpp(SEXP fSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP rSEXP, SEXP sSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(fourierin_fourierin_cx_2d_cpp_try(fSEXP, aSEXP, bSEXP, cSEXP, dSEXP, rSEXP, sSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// fft_rcpp
Rcpp::ComplexVector fft_rcpp(const Rcpp::NumericVector& real, const Rcpp::NumericVector& imag);
static SEXP fourierin_fft_rcpp_try(SEXP realSEXP, SEXP imagSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type real(realSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type imag(imagSEXP);
    __result = Rcpp::wrap(fft_rcpp(real, imag));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP fourierin_fft_rcpp(SEXP realSEXP, SEXP imagSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(fourierin_fft_rcpp_try(realSEXP, imagSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// fft_rcpp_2
Rcpp::ComplexVector fft_rcpp_2(const Rcpp::ComplexVector& v);
static SEXP fourierin_fft_rcpp_2_try(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const Rcpp::ComplexVector& >::type v(vSEXP);
    __result = Rcpp::wrap(fft_rcpp_2(v));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP fourierin_fft_rcpp_2(SEXP vSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(fourierin_fft_rcpp_2_try(vSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// fft_rcpp_3
Rcpp::ComplexVector fft_rcpp_3(const Rcpp::ComplexVector& v);
static SEXP fourierin_fft_rcpp_3_try(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const Rcpp::ComplexVector& >::type v(vSEXP);
    __result = Rcpp::wrap(fft_rcpp_3(v));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP fourierin_fft_rcpp_3(SEXP vSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(fourierin_fft_rcpp_3_try(vSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// fft_rcpp_4
arma::cx_vec fft_rcpp_4(const arma::cx_vec& v);
static SEXP fourierin_fft_rcpp_4_try(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type v(vSEXP);
    __result = Rcpp::wrap(fft_rcpp_4(v));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP fourierin_fft_rcpp_4(SEXP vSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(fourierin_fft_rcpp_4_try(vSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}

// validate (ensure exported C++ functions exist before calling them)
static int fourierin_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::cx_vec(*fourierin_1d_cpp)(const arma::vec&,double,double,double,double,double,double)");
        signatures.insert("arma::cx_vec(*fourierin_1d_nonregular_cpp)(const arma::vec&,double,double,const arma::vec&,int,double,double)");
        signatures.insert("arma::cx_mat(*fourierin_2d_cpp)(const arma::mat&,const arma::vec&,const arma::vec&,const arma::vec&,const arma::vec&,double,double)");
        signatures.insert("arma::cx_vec(*fourierin_cx_1d_cpp)(const arma::cx_vec&,double,double,double,double,double,double)");
        signatures.insert("arma::cx_vec(*fourierin_cx_1d_nonregular_cpp)(const arma::cx_vec&,double,double,const arma::vec&,int,double,double)");
        signatures.insert("arma::cx_mat(*fourierin_cx_2d_cpp)(const arma::cx_mat&,const arma::vec&,const arma::vec&,const arma::vec&,const arma::vec&,double,double)");
        signatures.insert("Rcpp::ComplexVector(*fft_rcpp)(const Rcpp::NumericVector&,const Rcpp::NumericVector&)");
        signatures.insert("Rcpp::ComplexVector(*fft_rcpp_2)(const Rcpp::ComplexVector&)");
        signatures.insert("Rcpp::ComplexVector(*fft_rcpp_3)(const Rcpp::ComplexVector&)");
        signatures.insert("arma::cx_vec(*fft_rcpp_4)(const arma::cx_vec&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP fourierin_RcppExport_registerCCallable() { 
    R_RegisterCCallable("fourierin", "fourierin_fourierin_1d_cpp", (DL_FUNC)fourierin_fourierin_1d_cpp_try);
    R_RegisterCCallable("fourierin", "fourierin_fourierin_1d_nonregular_cpp", (DL_FUNC)fourierin_fourierin_1d_nonregular_cpp_try);
    R_RegisterCCallable("fourierin", "fourierin_fourierin_2d_cpp", (DL_FUNC)fourierin_fourierin_2d_cpp_try);
    R_RegisterCCallable("fourierin", "fourierin_fourierin_cx_1d_cpp", (DL_FUNC)fourierin_fourierin_cx_1d_cpp_try);
    R_RegisterCCallable("fourierin", "fourierin_fourierin_cx_1d_nonregular_cpp", (DL_FUNC)fourierin_fourierin_cx_1d_nonregular_cpp_try);
    R_RegisterCCallable("fourierin", "fourierin_fourierin_cx_2d_cpp", (DL_FUNC)fourierin_fourierin_cx_2d_cpp_try);
    R_RegisterCCallable("fourierin", "fourierin_fft_rcpp", (DL_FUNC)fourierin_fft_rcpp_try);
    R_RegisterCCallable("fourierin", "fourierin_fft_rcpp_2", (DL_FUNC)fourierin_fft_rcpp_2_try);
    R_RegisterCCallable("fourierin", "fourierin_fft_rcpp_3", (DL_FUNC)fourierin_fft_rcpp_3_try);
    R_RegisterCCallable("fourierin", "fourierin_fft_rcpp_4", (DL_FUNC)fourierin_fft_rcpp_4_try);
    R_RegisterCCallable("fourierin", "fourierin_RcppExport_validate", (DL_FUNC)fourierin_RcppExport_validate);
    return R_NilValue;
}
