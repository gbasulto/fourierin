
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fourierin

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/gbasulto/fourierin?branch=master&svg=true)](https://ci.appveyor.com/project/gbasulto/fourierin)
[![Travis build
status](https://travis-ci.org/gbasulto/fourierin.svg?branch=master)](https://travis-ci.org/gbasulto/fourierin)
<!-- badges: end -->

This is a package in `R` to numerically calculate Fourier-type integrals
of univariate and bivariate functions with compact support and
simultaneously evaluated at several points. If the evaluation grid is
equally spaced, a Fast Fourier Transform method is used to speed up
computations.

See detailed documentation on the vignette (type
`browseVignettes("fourierin")` and then click “HTML”).

## Installation

You can install the released version of fourierin from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fourierin")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gbasulto/fourierin")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fourierin)
## basic example code
```
