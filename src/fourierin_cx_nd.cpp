#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

// using namespace Rcpp;

/*

  Fast Fourier Transform for Rcpp complex vectors

  This function performs the fft by transforming it to a RcppArmadillo
  vector ann then use the function in this library.

*/

//' @export
//[[Rcpp::export]]
Rcpp::ComplexVector fft_rcpp(Rcpp::NumericVector real, Rcpp::NumericVector imag)
{
  Rprintf("v = %f\n", v[0].i);
  // // Reuses memory and avoids extra copy
  // arma::vec re(), im;
  // arma::cx_vec v_arma(v.begin(), v.size());

  // //  arma::fft();
  // v_arma.print("v = ");
    
  return v;
}
