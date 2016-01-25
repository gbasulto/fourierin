#' Compute Fourier integrals
#'
#' It computes Fourier integrals for functions of one and two
#' variables.
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

