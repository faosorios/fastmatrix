## ID: mahalanobis.R, last updated 2020-10-23, F.Osorio

Mahalanobis <- function(x, center, cov, inverted = FALSE)
{ ## computes the squared Mahalanobis distance
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  if (is.null(center))
    stop("'center' must be provided")
  if (isFALSE(center))
    center <- double(p) # center is zeroed
  if (!is.vector(center))
    stop("'center' must be a vector")
  if (length(center) != p)
    stop("'center' has incorrect length")

  if (is.null(cov))
    stop("'cov' matrix must be provided")
  if (!is.matrix(cov))
    stop("'cov' must be a matrix")
  if (all(dim(cov) != c(p,p)))
    stop("'cov' has incorrect dimensions")

  storage.mode(x) <- "double"
  storage.mode(cov) <- "double"

  z <- .C("mahal_distances",
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          center = as.double(center),
          cov = cov,
          inverted = as.integer(inverted),
          distances = double(n))$distances
  z
}
