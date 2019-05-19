context("univariate")

test_that("fourierin works for 1d for real function", {
  ## Function to be used in the integrand
  myfnc <- function(t) exp(-t^2/2)

  ## Compute integral
  out <- fourierin(f = myfnc, lower_int = -5, upper_int = 5,
                   lower_eval= -3, upper_eval = 3, const_adj = -1,
                   freq_adj = -1, resolution = 64)
  expect_is(out, "list")
})

test_that("fourierin works for 1d for complex function", {
  ## Function to be used in the integrand
  myfnc <- function(t) exp(-1i*t^2/2)

  ## Compute integral
  out <- fourierin(f = myfnc, lower_int = -5, upper_int = 5,
                   lower_eval= -3, upper_eval = 3, const_adj = -1,
                   freq_adj = -1, resolution = 64)
  expect_is(out, "list")
})


test_that(
  "error if evaluation grid nor evaluation bounds are provided", {
  ## Function to be used in the integrand
  myfnc <- function(t) exp(-t^2/2)

  ## Compute integral
  expect_error(
    fourierin(f = myfnc, lower_int = -5, upper_int = 5,
              const_adj = -1,
              freq_adj = -1, resolution = 64))
})


test_that(
  "evaluate on a given set of points (real function)", {
    myfnc <- function(t) exp(-t^2/2)
    out <- fourierin(f = myfnc, lower_int = -5, upper_int = 5,
                     const_adj = -1, freq_adj = -1, resolution = 64,
                     eval_grid = 0:1)
    expect_type(out, "complex")
  }
)

test_that(
  "evaluate on a given set of points (complex function)", {
    myfnc <- function(t) exp(-1i*t^2/2)
    out <- fourierin(f = myfnc, lower_int = -5, upper_int = 5,
                     const_adj = -1, freq_adj = -1, resolution = 64,
                     eval_grid = 0:1)
    expect_type(out, "complex")
  }
)


test_that(
  "evaluate f outside fourierin", {
    t <- seq(-5, 5 - 10/64, length.out = 64) + 10/64/2
    f <- exp(-t^2/2)
    out <- fourierin(f = f, lower_int = -5, upper_int = 5,
                     const_adj = -1, freq_adj = -1,
                     eval_grid = 0:1)

    expect_type(out, "complex")
  }
)









