
context("bivariate")

test_that("fourierin works for 2d", {
  ##-Parameters of bivariate normal distribution
  mu <- c(-1, 1)
  sig <- matrix(c(3, -1, -1, 2), 2, 2)
  
  ##-Multivariate normal density
  ##-x is n x d
  f <- function(x) {
    ##-Auxiliar values
    d <- ncol(x)
    z <- sweep(x, 2, mu, "-")
    ##-Get numerator and denominator of normal density
    num <- exp(-0.5*rowSums(z * (z %*% solve(sig))))
    denom <- sqrt((2*pi)^d*det(sig))
    return(num/denom)
  }
  
  
  ##-Approximate cf using Fourier integrals
  eval <- fourierin(f, lower_int = c(-8, -6), upper_int = c(6, 8),
                    lower_eval = c(-4, -4), upper_eval = c(4, 4),
                    const_adj = 1, freq_adj =  1,
                    resolution = c(128, 128))
  expect_is(eval, "list")
})


