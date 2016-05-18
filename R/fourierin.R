#' Univariate Fourier integrals
#'
#' It computes Fourier integrals of functions of one and two
#' variables on a regular grid.
#'
#' @param f function or a vector of size m. If a function is provided,
#'     it must be able to be evaluated in vectors. If a vector of
#'     values is provided, such evaluations must have been obtained on
#'     a regular grid and the Fourier integral is faster is m is a
#'     power of 2.
#' @param a Lower integration limit.
#' @param b Upper integration limit.
#' @param c Lower evaluation limit.
#' @param d Upper evaluation limit.
#' @param r Factor related to adjust definition of Fourier
#'     transform. It is usually equal to 0, -1 or 1.
#' @param s Constant to adjust the exponent on the definition of the
#'     Fourier transform. It is usually equal to 1, -1, 2pi or -2pi.
#' @param resol An integer (faster if power of two) determining the
#'     resolution of the evaluation grid. Not required if f is a
#'     vector.
#' @param w An optional vector where the Fourier integral will be
#'     evaluated. If provided, the FFT will NOT be used.
#' @param use_fft Logical value specifying whether the FFT will be
#'     used to compute the Fourier integral.
#' @return A list with the elements \item{w}{A vector of size m where
#'     the integral was computed.}  \item{values}{A complex vector of
#'     size m with the values of the integral}
#'
#' @example
#' examples/ex_fourierin_1d.R
#' @export
fourierin_1d <- function(f, a, b, c = NULL, d = NULL,
                         r, s, resol = NULL, w = NULL,
                         use_fft = TRUE) {
    ## If function values are provided, then the resolution
    ## is the length of the vector of values.
    if (!is.function(f)) resol <- length(f)

    ## Increment in the frequency domain.
    gam <- (d - c)/resol

    ## Freq. dom. vector. If w is provided, FFT will NOT be used.
    if (is.null(w)) {
        w <- seq(c, d - gam, length.out = resol)
        if (is.null(c) | is.null(d)) {
            stop("c and d must be provided.")
        }
    } else {
        use_fft <- FALSE
    }

    ## If f is the function, it needs to be evaluated in
    ## the time domain values.
    if (is.function(f)) {

        del <- (b - a)/resol # Increment in the time domain.
        t <- seq(a + del/2, b - del/2,
                 length.out = resol)    # Freq. dom. vector.
        f_t <- f(t)                     # Function values
        ## Rutinary check
        if(is.null(f_t)) stop("Function f is null.")
    } else {
        f_t <- f
    }

    if (!use_fft) {
        out <- switch(is.complex(f_t) + 1,
                      fourierin_1d_nonregular_cpp(f_t, a, b, w,
                                                  resol, r, s),
                      fourierin_cx_1d_nonregular_cpp(f_t,
					 a, b, w, resol, r, s))
    } else {
        out <- switch(is.complex(f_t) + 1,
                      fourierin_1d_cpp(f_t, a, b, c, d, r, s),
                      fourierin_cx_1d_cpp(f_t, a, b, c, d, r, s))
    }

    return(list(w = w,                  # Return list.
                values = out))
}

#' Bivariate Fourier integrals
#'
#' It computes Fourier integrals for functions of one and two
#' variables.
#'
#' @param f function or a matrix of size m1 x m2. If a function is
#'     provided, it must be able to be evaluated in a matrix of two
#'     columns. If a matrix of values is provided instead, such
#'     evaluations must have been obtained on a regular grid matrix
#'     and the Fourier integral is faster is m1 and m2 are powers of
#'     2.
#' @param a Lower integration limits.
#' @param b Upper integration limits.
#' @param c Lower evaluation limits.
#' @param d Upper evaluation limits.
#' @param r Factor related to adjust definition of Fourier
#'     transform. It is usually 0, -1 or 1.
#' @param s Constant to adjust the exponent on the definition of the
#'     Fourier transform. It is usually 1, -1, 2pi or -2pi.
#' @param resol A vector of two integers (faster if entries are powers
#'     of two) determining the resolution of the evaluation grid. Not
#'     required if f is a vector.
#'
#' @return A list with three elements
#' \item{w1}{Evaluation grid for first entry}
#' \item{w2}{Evaluation grid for second entry}
#' \item{values}{m1 x m2 matrix of complex numbers, correspoding to the
#'       evaluations of the integral}
#'
#' @example
#' examples/ex_fourierin_2d.R
#' @export
fourierin_2d <- function(f, a, b, c, d, r, s, resol = NULL,
                         w = NULL, use_fft = TRUE){
    ## If function values are provided, then the resolution is the
    ## length of the vector of values.
    if(!is.function(f)) resol <- dim(f)

    ## Increment in the frequency domain.
    gam <- (d - c)/resol

    ## Freq. dom. vectors.
    if (is.null(w)){
        w1 <- seq(c[1], d[1] - gam[1], length.out = resol[1])
        w2 <- seq(c[2], d[2] - gam[2], length.out = resol[2])
    }

    ## If f is the function, it needs to be evaluated in the time
    ## domain values.
    if(is.function(f)){
        del <- (b - a)/resol # Increment in the time domain.
        t1 <- seq(a[1] + del[1]/2, b[1] - del[1]/2,
                  length.out = resol[1]) # Freq. dom. vector.
        t2 <- seq(a[2] + del[2]/2, b[2] - del[2]/2,
                  length.out = resol[1]) # Freq. dom. vector.
        t <- as.matrix(expand.grid(t1, t2))
        f_vals <- matrix(f(t), resol[1], resol[2])

                                        # Rutinary check
        if(is.null(f_vals)) stop("Function f is null.")

    } else{
        f_vals <- f
    }

    if (!use_fft) {
        w_temp <- switch(is.null(w) + 1,
                         w,
                         as.matrix(expand.grid(w1, w2))
                         )

        out <- switch(is.complex(f_vals) + 1,
                      fourierin_2d_nonregular_cpp(f_vals, a, b,
                                                  w_temp,
                                                  resol, r, s),
                      fourierin_cx_2d_nonregular_cpp(f_vals, a, b,
                                                     w_temp, resol,
                                                     r, s))
        ## If no grid was provided, the values of the Fourier integral
        ## are put into a matrix.
        if(is.null(w)) out <- matrix(out, resol[1], resol[2])
    } else {
        out <- switch(is.complex(f_vals) + 1,
                      fourierin_2d_cpp(f_vals, a, b, c, d, r, s),
                      fourierin_cx_2d_cpp(f_vals, a, b, c, d, r, s))
    }


    ## Return grid if it was not provided.
    if (is.null(w)) {
      return(list(w1 = w1,
                  w2 = w2,
                  values = out))

    } else {
      return(out)
    }
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
#' @param resol An integer (faster if power of two) determining the
#'     resolution of the evaluation grid. Not required if f is a vector.
#'
#' @return A list with the elements n-dimensional array and n vectors
#'     with their corresponding resolution. Specifically,
#' \item{values}{A n-dimensional (resol_1 x resol_2 x ... x resol_n)
#'       complex array with the values.}
#' \item{w1}{A vector of size resol_1}
#' \item{...}{ }
#' \item{wn}{A vector of size resol_n}
#'
#' @example
#' examples/ex_fourierin.R
#' @export
fourierin <- function(f, a, b, c, d, r, s, resol = NULL,
                      w = NULL, use_fft = TRUE){

  n <- length(a)                      # Get dimension of function
  # from lower integration
  # limit.
  switch(n,
         return(fourierin_1d(f, a, b, c, d, r, s, resol, w, use_fft)),
         return(fourierin_2d(f, a, b, c, d, r, s, resol, w, use_fft))
  ) # End switch
}

#' @useDynLib fourierin
#' @importFrom Rcpp sourceCpp
NULL

