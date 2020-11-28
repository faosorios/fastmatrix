## ID: mardia.R, last updated 2020-11-22, F.Osorio

kurtosis <- function(x)
{ ## Mardia's multivariate kurtosis coefficient
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  z <- cov.weighted(x, weights = rep(1, n))
  storage.mode(x) <- "double"
  storage.mode(z$cov) <- "double"

  z <- .C("skewness_and_kurtosis",
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          center = as.double(z$mean),
          cov = z$cov,
          stats = double(2),
          task = as.integer(0))$stats[2]
  attr(z, 'excess') <- z - p * (p + 2)
  z
}

skewness <- function(x)
{ ## Mardia's multivariate skewness coefficient
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  z <- cov.weighted(x, weights = rep(1, n))
  storage.mode(x) <- "double"
  storage.mode(z$cov) <- "double"

  z <- .C("skewness_and_kurtosis",
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          center = as.double(z$mean),
          cov = z$cov,
          stats = double(2),
          task = as.integer(1))$stats[1]
  z
}
