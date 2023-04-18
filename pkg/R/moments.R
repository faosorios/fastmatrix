## ID: moments.R, last updated 2023-04-17, F.Osorio

moments <- function(x)
{ ## calculates up to fourth central moments
  if (!is.vector(x))
    stop("argument 'x' must be a vector.")
  if (!is.numeric(x))
    stop("argument 'x' is not a numeric vector.")

  # removing NA's
  ok <- !is.na(x)
  x <- x[ok]

  n <- length(x)
  z <- .C("central_moments",
          x = as.double(x),
          n = as.integer(n),
          mean = double(1),
          m2 = double(1),
          m3 = double(1),
          m4 = double(1))[c("mean","m2","m3","m4")]
  z <- list(second = z$m2, third = z$m3, fourth = z$m4, skewness = z$m3 / z$m2^(1.5),
            kurtosis = z$m4 / z$m2^2)
  z
}

geomean <- function(x)
{ ## geometric mean using a compensated product scheme
  if (!is.vector(x))
    stop("argument 'x' must be a vector.")
  if (!is.numeric(x))
    stop("argument 'x' is not a numeric vector.")

  # removing NA's
  ok <- !is.na(x)
  x <- x[ok]

  if (any(x <= 0)) {
    warning("Non-positive values in 'x'.")
    return(NA)
  }

  n <- length(x)
  y <- x^(1 / n)
  z <- .C("geometric_mean",
          y = as.double(y),
          n = as.integer(n),
          mean = double(1))$mean
  z
}
