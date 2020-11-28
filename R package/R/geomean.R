## ID: geomean.R, last updated 2020-11-27, F.Osorio

geomean <- function(x)
{ ## geometric mean using a compensated product scheme
  if (!is.vector(x))
    stop("argument 'x' must be a vector.")
  if (!is.numeric(x))
    stop("argument 'x' is not a numeric vector.")

  if (any(x <= 0)) {
    warning("Non-positive values in 'x'.")
    return(NA)
  }

  n <- length(x)
  z <- .C("geometric_mean",
          x = as.double(x),
          n = as.integer(n),
          mean = double(1))$mean
  z
}
