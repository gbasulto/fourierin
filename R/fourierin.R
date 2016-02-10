#' Compute Univariate Fourier integrals
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
#' examples/ex_fourierin_1d.R
#' examples/ex_fourierin_2d.R
#' @export
fourierin_1d <- function(f, a, b, c, d, r, s, resol = NULL){
    ## If function values are provided, then the resolution
    ## is the length of the vector of values.
    if(!is.function(f)) resol <- length(f)

    ## Increment in the frequency domain.
    gam <- (d - c)/resol

    ## Freq. dom. vector.
    w <- seq(c, d - gam, length.out = resol)

    ## If f is the function, it needs to be evaluated in
    ## the time domain values.
    if(is.function(f)){
        del <- (b - a)/resol # Increment in the time
                                        # domain.
        t <- seq(a + del/2, b - del/2,
                 length.out = resol)    # Freq. dom. vector.
        out <- fourierin_1d_cpp(f(t), a, b, c, d, r, s)
    } else{
        out <- fourierin_1d_cpp(f, a, b, c, d, r, s)
    }

    return(list(w = w,                  # Return list.
                values = out))
}

#' Compute Univariate Fourier integrals
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
#' examples/ex_fourierin_1d.R
#' @export
fourierin_2d <- function(f, a, b, c, d, r, s, resol = NULL){
    ## If function values are provided, then the resolution is the
    ## length of the vector of values.
    if(!is.function(f)) resol <- dim

    ## Increment in the frequency domain.
    gam <- (d - c)/resol

    ## Freq. dom. vectors.
    w1 <- seq(c[1], d[1] - gam[1], length.out = resol[1])
    w2 <- seq(c[2], d[2] - gam[2], length.out = resol[2])

    ## If f is the function, it needs to be evaluated in the time
    ## domain values.
    if(is.function(f)){
        del <- (b - a)/resol # Increment in the time
                                        # domain.
        t1 <- seq(a[1] + del[1]/2, b[1] - del[1]/2,
                  length.out = resol[1]) # Freq. dom. vector.
        t2 <- seq(a[2] + del[2]/2, b[2] - del[2]/2,
                  length.out = resol[1]) # Freq. dom. vector.
        t <- as.matrix(expand.grid(t1, t2))
        f_vals <- matrix(f(t), resol[1], resol[2])
        out <- fourierin_2d_cpp(f_vals, a, b, c, d, r, s)
#         out <- fourierin_2d_cpp(outer(t1, t2, f), a, b, c, d, r, s)
    } else{
        out <- fourierin_2d_cpp(f, a, b, c, d, r, s)
    }

    return(list(w1 = w1,
                w2 = w2,
                values = out))
}



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
#' examples/ex_fourierin_1d.R
#' @export
fourierin <- function(f, a, b, c, d, r, s, resol = NULL){

    n <- length(a)                      # Get dimension of function
                                        # from lower integration
                                        # limit.
    switch(n,
           return(fourierin_1d(f, a, b, c, d, r, s, resol)),
           return(fourierin_2d(f, a, b, c, d, r, s, resol))
           ) # End switch
}

#' @useDynLib fourierin
#' @importFrom Rcpp sourceCpp
NULL

