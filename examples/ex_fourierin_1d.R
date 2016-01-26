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
                 a = -15, b = 15, c = 0, d = 8,
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

