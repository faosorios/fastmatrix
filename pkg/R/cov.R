## ID: cov.weighted.R, last updated 2020-09-23, F.Osorio

cov.weighted <- function(x, weights = rep(1, nrow(x)))
{ ## computes weighted mean and covariance using an online algorithm
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)
  if (!missing(weights)) {
    if (length(weights) != n)
      stop("length of 'weights' must equal the number of rows in 'x'")
    if (any(weights < 0) || (sum(weights) == 0))
	    stop("weights must be non-negative and not all zero")
  }

  storage.mode(x) <- "double"
  y <- .C("cov_weighted",
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          weights = as.double(weights),
          mean = double(p),
          cov  = double(p * p))[c("mean","cov")]
  y$cov <- matrix(y$cov, nrow = p, ncol = p)
  y
}

cov.MSSD <- function(x)
{ ## covariance estimation using the Mean Square Successive Method (MSSD)
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  storage.mode(x) <- "double"
  y <- .C("cov_MSSD",
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          mean = double(p),
          cov  = double(p * p))[c("mean","cov")]
  y$cov <- matrix(y$cov, nrow = p, ncol = p)
  y
}
