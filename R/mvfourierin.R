
## -------------------------------------------------------------------

## The next functions are only a test.
z <- function (m, delta) {
    
    ## Argument
    idx <- 0:(m - 1)
    arg1 <- -delta*idx^2
    arg2 <- -delta*(idx - m)^2
    z <- complex(argument = c(arg1, arg2))
    
    fft(z)
}
