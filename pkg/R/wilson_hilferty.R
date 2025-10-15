## ID: wilson_hilferty.R, last updated 2024-01-03, F.Osorio

wilson.hilferty <- function(x, shape, rate = 1)
{ ## Wilson-Hilferty transformation for gamma deviates
  if (!is.vector(x))
    stop("'x' must be a vector")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  if (any(x < 0)) 
    stop("'x' must be non-negative")
  n <- length(x)

  if (is.null(shape))
    stop("'shape' must be provided")
  if (shape < 0) 
    stop("'shape' must be non-negative")
  if (rate < 0) 
    stop("'rate' must be non-negative")

  z <- .C("wilson_hilferty_gamma",
          x = as.double(x),
          n = as.integer(n),
          shape = as.double(shape),
          rate = as.double(rate),
          z = double(n))$z
  z
}

WH.normal <- function(x)
{ ## Wilson-Hilferty transformation for multivariate normal deviates
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  z <- cov.weighted(x, weights = rep(1, n))
  D2 <- Mahalanobis(x, center = z$mean, cov = z$cov, inverted = FALSE)
  z <- .C("wilson_hilferty_chisq",
          distances = as.double(D2),
          n = as.integer(n),
          p = as.integer(p),
          z = double(n))$z
  z
}
