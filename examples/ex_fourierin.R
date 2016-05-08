## -- Example 0 ------------------------------------------------------
##                       Speed test 
## -------------------------------------------------------------------
library(fourierin)

## Test speed at several resolutions
resolution <- 2^(0:8)
## Compute the time for every resolution
times <- 
    t(apply(matrix(resolution), 1,
            function(resol) {
                out <- microbenchmark::microbenchmark(
                    fourierin_1d(f = function(t) exp(-t^2/2),
                                 -5, 5, -3, 3, -1, -1, resol),
                    fourierin_1d(f = function(t) exp(-t^2/2),
                                 -5, 5, -3, 3, -1, -1, resol,
                                 use_fft = FALSE),
                    times = 5
                )
                out <- with(out, tapply(time/1000, expr, mean))
                out <- c(out[[1]], out[[2]])
            }
            )
      )
## Convert to a dataframe with appropriate names
summ <- data.frame(resolution,
                   FFT = times[, 1], no_FFT = times[, 2])
## ... And plot it
with(summ, {
    plot(range(resolution), range(c(FFT, no_FFT)), type = "n",
         xlab = "resolution", ylab = "time (microseconds)")
    lines(resolution, FFT, type = "b", col = "cyan")
    lines(resolution, no_FFT, type = "b", col = "magenta")
    legend("topleft", legend = c("FFT", "No FFT"),
           col = c("cyan", "magenta"), lty = 1, pch = 1)
})


## -- Example 1 ------------------------------------------------------
## -- Recovering std. normal from its characteristic function --------
library(fourierin)

                                        # Compute integral
out <- fourierin(f = function(t) exp(-t^2/2),
                 a = -5, b = 5, c = -3, d = 3,
                 r = -1, s = -1, resol = 64)
grid <- out$w                           # Extract grid and values
values <- Re(out$values)

plot(grid, values, type = "l", col = 3)
lines(grid, dnorm(grid), col = 4)

## -- Example 2 ------------------------------------------------------
## -- Computing characteristic function of a gamma r. v. ------------

library(fourierin)
                                        # Compute integral
shape <- 5
rate <- 3
out <- fourierin(f = function(t) dgamma(t, shape, rate),
                 a = 0, b = 6, c = -4, d = 4,
                 r = 1, s = 1, resol = 64)
grid <- out$w                           # Extract grid
re_values <- Re(out$values)             # Real values
im_values <- Im(out$values)             # Imag values

                                        # Now compute the real and
                                        # imaginary true values of the
                                        # characteric function.
true_cf <- function(t, shape, rate) (1 - 1i*t/rate)^-shape
true_re <- Re(true_cf(grid, shape, rate))
true_im <- Im(true_cf(grid, shape, rate))

                                        # Compare them. We can see a
                                        # slight discrepancy on the
                                        # tails, but that is fixed
                                        # when resulution is
                                        # increased.
plot(grid, re_values, type = "l", col = 3)
lines(grid, true_re, col = 4)

                                        # Same here
plot(grid, im_values, type = "l", col = 3)
lines(grid, true_im, col = 4)

## -- Example 3 ------------------------------------------------------
## -- Recovering std. normal from its characteristic function --------
library(fourierin)

## Parameters of bivariate normal distribution
mu <- c(-1, 1)
sig <- matrix(c(3, -1, -1, 2), 2, 2)

## Multivariate normal density
## x is n x d
f <- function(x) {
    ## Auxiliar values
    d <- ncol(x)
    z <- sweep(x, 2, mu, "-")

    ## Get numerator and denominator of normal density
    num <- exp(-0.5*rowSums(z * (z %*% solve(sig))))
    denom <- sqrt((2*pi)^d*det(sig))

    return(num/denom)
}

## Characteristic function
## s is n x d
phi <- function(s) {
    complex(modulus = exp(- 0.5*rowSums(s*(s %*% sig))),
            argument = s %*% mu)
}

## Approximate cf using Fourier integrals
eval <- fourierin(f, a = c(-8, -6), b = c(6, 8),
                    c = c(-4, -4), d = c(4, 4),
                    r = 1, s = 1, resol = c(128, 128))
t1 <- eval$w1
t2 <- eval$w2
t <- as.matrix(expand.grid(t1 = t1, t2 = t2))
approx <- eval$values
true <- matrix(phi(t), 128, 128)        # Compute true values

## This is a section of the characteristic functions
i <- 65
plot(t2, Re(approx[i, ]), type = "l", col = 2,
     ylab = "",
     xlab = expression(t[2]),
     main = expression(paste("Real part section at ",
                             t[1], "= 0")))
lines(t2, Re(true[i, ]), col = 3)
legend("topleft", legend = c("true", "approximation"),
       col = 3:2, lwd = 1)

## Another section, now of the imaginary part
plot(t1, Im(approx[, i]), type = "l", col = 2,
     ylab = "",
     xlab = expression(t[1]),
     main = expression(paste("Imaginary part section at ",
                             t[2], "= 0")))
lines(t1, Im(true[, i]), col = 3)
legend("topleft", legend = c("true", "approximation"),
       col = 3:2, lwd = 1)


