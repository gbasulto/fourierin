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


