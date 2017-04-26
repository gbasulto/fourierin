#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP fourierin_fourierin_1d_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fourierin_fourierin_1d_nonregular_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fourierin_fourierin_2d_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fourierin_fourierin_2d_nonregular_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fourierin_fourierin_cx_1d_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fourierin_fourierin_cx_1d_nonregular_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fourierin_fourierin_cx_2d_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fourierin_fourierin_cx_2d_nonregular_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fourierin_RcppExport_registerCCallable();

static const R_CallMethodDef CallEntries[] = {
  {"fourierin_fourierin_1d_cpp",               (DL_FUNC) &fourierin_fourierin_1d_cpp,               7},
  {"fourierin_fourierin_1d_nonregular_cpp",    (DL_FUNC) &fourierin_fourierin_1d_nonregular_cpp,    7},
  {"fourierin_fourierin_2d_cpp",               (DL_FUNC) &fourierin_fourierin_2d_cpp,               7},
  {"fourierin_fourierin_2d_nonregular_cpp",    (DL_FUNC) &fourierin_fourierin_2d_nonregular_cpp,    7},
  {"fourierin_fourierin_cx_1d_cpp",            (DL_FUNC) &fourierin_fourierin_cx_1d_cpp,            7},
  {"fourierin_fourierin_cx_1d_nonregular_cpp", (DL_FUNC) &fourierin_fourierin_cx_1d_nonregular_cpp, 7},
  {"fourierin_fourierin_cx_2d_cpp",            (DL_FUNC) &fourierin_fourierin_cx_2d_cpp,            7},
  {"fourierin_fourierin_cx_2d_nonregular_cpp", (DL_FUNC) &fourierin_fourierin_cx_2d_nonregular_cpp, 7},
  {"fourierin_RcppExport_registerCCallable",   (DL_FUNC) &fourierin_RcppExport_registerCCallable,   0},
  {NULL, NULL, 0}
};

void R_init_fourierin(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
