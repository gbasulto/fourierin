#' Compute Fourier integrals
#'
#' It computes Fourier integrals for functions of one and two
#' variables.
#'
#' @param f A function which can be evaluated in matrices of n columns
#'     (n = 1 or n = 2). Or a matrix of n columns with f already
#'     evaluated.
#' @param a Lower integration limit(s).
#' @param b Upper integration limit(s).
#' @param c Lower evaluation limit(s).
#' @param d Upper evaluation limit(s).
#' @param r Factor related to adjust definition of Fourier
#'     transform. It is usually 0, -1 or 1.
#' @param s Constant to adjust the exponent on the definition of the
#'     Fourier transform. It is usually 1, -1, 2pi or -2pi.
#'
#' @return A list with the elements
#' \item{w}{ddd}
#' \item{values}{dddd}
#'
#' @example
#' examples/fourierin_1d.R
#' @export
fourierin <- function(f, a, b, c, d, r, s, resol){
                                        # If function values are
                                        # provided, then the
                                        # resolution is the length of
                                        # the vector of values.
    if(!is.function(f)) resol <- length(f)

    gam <- (d - c)/resol                # Increment in the frequency
                                        # domain.

    w <- seq(c, d - gam, length.out = resol) # Freq. dom. vector.

                                        # If f is the function, it
                                        # needs to be evaluated in the
                                        # time domain values.
    if(is.function(f)){
        del <- (b - a)/resol            #Increment in the time
                                        # domain.
        t <- seq(a, b - del, length.out = resol) # Freq. dom. vector.
        out <- fourierin_1d(f(t), a, b, c, d, r, s)
    } else{
        out <- fourierin_1d(f, a, b, c, d, r, s)
    }

    return(list(w = w,                  # Return list.
                values = out))
}

#' @useDynLib fourierin
#' @importFrom Rcpp sourceCpp
NULL

