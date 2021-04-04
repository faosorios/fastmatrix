## ID: whitening.R, last updated 2021-04-04, F.Osorio

whitening <- function(x, Scatter = NULL)
{ # Cholesky whitening transformation
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  storage.mode(x) <- "double"
  n <- nrow(x)
  p <- ncol(x)

  if (is.null(Scatter))
    Scatter <- cov.weighted(x, weights = rep(n / (n - 1), n))$cov

  if (!is.matrix(Scatter))
    stop("'Scatter' must be a matrix")
  if (nrow(Scatter) != ncol(Scatter))
    stop("'Scatter' must be a square matrix")
  if (p != ncol(Scatter))
    stop("'x' and 'Scatter' are not conformable")
  storage.mode(Scatter) <- "double"

  z <- .C("whitening_chol",
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          Scatter = Scatter)$x
  z
}
