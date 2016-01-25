
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

arma::cx_vec fourierin_1d(arma::vec f, double a,
			  double b, double c, double d, double r)
{
  int m = f.n_rows;
  arma::cx_vec out(m), y(2*m), z(2*m), aux(2*m);
  arma::vec J1(m), J2(m), t(m), w(m), arg(m);
  double bet, gam, del, cnst;

  bet = (b - a)/m;		// Real numbers
  gam = (d - c)/m;
  del = bet*gam/2;
  J1 = arma::linspace<arma::vec>(0, m-1, m); // mx1 vectors
  J2 = arma::linspace<arma::vec>(m, 2*m-1, m);
  t = a + bet*J1;
  w = c + gam*J1;
  y.zeros();		// (2m) x 1 vector

  // We will first compute the argument and then create the complex
  // vector.
  arg = J1 % (bet*c + del*J1);	// Fill y
  y.rows(0, m - 1) = cx_vec(f % cos(arg), f % sin(arg));

  arg = -del*pow(J1, 2);        // Fill first half of z
  z.rows(0, m - 1) = cx_vec(cos(arg), sin(arg));

  arg = -del*pow(J2 - 2*m, 2);  // Fill 2nd. half of z
  z.rows(m, 2*m - 1) = cx_vec(cos(arg), sin(arg));

  aux = ifft(fft(y) % fft(z));
  out = aux.rows(0, m - 1);
  cnst = bet*pow(2*datum::pi, -(1 - r)/2);
  arg = (a*w + del*pow(J1, 2));
  out = cnst*(out % cx_vec(cos(arg), sin(arg)));

  return out;
}

//' Computes Fourier integral of univariate functions
//'
//' This function computes univariate and bivariate continuous
//' Fourier tranform based on the paper by Inverarity (2002):
//' "Fast computation of multidimensional Fourier integrals".
//' It is the formula (4.1) on the paper.
//'
//' @param f Values of the function.
//' @param m Resolution of the integral.
//' @param a \eqn{n x 1} vector. Lower integration limit.
//' @param b \eqn{n x 1} vector. Upper integration limit.
//' @param c \eqn{n x 1} vector. Lower limit of w.
//' @param d \eqn{n x 1} vector. Upper limit of w.
//' @param r Power in (4.1).
//' @param s Scale constant in (4.1).
//' @return Continuous Fourier transform values at w.
//'
//' @export
// [[Rcpp::export]]
arma::cx_vec fourierin_1d(arma::vec f, double a,
			  double b, double c, double d,
			  double r, double s)
{

  int m = f.n_rows;
  arma::cx_vec out(m);

  // fourierin_1d without s argument is meant for s = 1. Thus we have
  // to make it valid for any s.
  out = pow(abs(s), 1/2)*fourierin_1d(f, a, b, s*c, s*d, r);

  return out;
}
