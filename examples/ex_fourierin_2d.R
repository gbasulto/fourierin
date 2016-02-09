## -- Example 3 ------------------------------------------------------
## -- Recovering std. normal from its characteristic function --------
library(fourierin)
library(MASS)                           # For biv. normal

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

# library(mvtnorm)
# library(ggplot2)
# library(dplyr)
# library(tidyr)
#
#
# x1 <- seq(-5, 5, len = 150) + mu[1]
# x2 <- seq(-5, 5, len = 150) + mu[2]
# expand.grid(x1 = x1, x2 = x2) %>%
#     tbl_df() %>%
#     mutate(y = f(cbind(x1, x2))) %>%
#     ggplot(aes(x1, x2)) +
#     geom_raster(aes(fill = y), interpolate = T)
#
#
# t1 <- seq(-2, 2, len = 150)
# t2 <- seq(-2, 2, len = 150)
# expand.grid(t1 = t1, t2 = t2) %>%
#     tbl_df() %>%
#     mutate(y = Im(phi(cbind(t1, t2)))) %>%
#     ggplot(aes(t1, t2)) +
#     geom_raster(aes(fill = y), interpolate = T)
#
# zz <- fourierin(f,
#                 a = c(-8, -6), b = c(6, 8),
#                 c = c(-2, -2), d = c(2, 2),
#                 r = 1, s = 1, resol = c(256, 256))
#
# expand.grid(t1 = zz$w1, t2 = zz$w2) %>%
#   tbl_df() %>%
#   mutate(y = Im(c(zz$values))) %>%
#   ggplot(aes(t1, t2)) +
#   geom_raster(aes(fill = y), interpolate = T)
#
#
# expand.grid(t1 = zz$w1, t2 = zz$w2) %>%
#   tbl_df() %>%
#   mutate(y = Re(c(zz$values))) %>%
#   filter(y > .773)


## expand.grid(x1 = x1, x2 = x2) %>%
##     tbl_df() %>%
##     mutate(y = dmvnorm(cbind(x1, x2), mu, sig)) %>%
##     ggplot(aes(x1, x2)) +
##     geom_raster(aes(fill = y), interpolate = T)


## microbenchmark(
##     dmvnorm(rbind(sig, sig), mu, sig),
##     f(rbind(sig, sig))
## )
