## ID: cholesky.R, last updated 2020-09-23, F.Osorio

chol2inverse <- function(x)
{ ## compute weighted mean and covariance using an online algorithm
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix")
  if (!is.upper.tri(x))
    stop("argument x is not an upper triangular matrix")

  storage.mode(x) <- "double"
  ldx <- p <- nrow(x)
  job <- 1 # 'x' is upper triangular

  y <- .C("chol_inverse",
          x = x,
          ldx = as.integer(ldx),
          p = as.integer(p),
          job = as.integer(job),
          info = as.integer(0))

  if (y$info > 0)
    stop("inverse of 'x' could not be computed")

  y <- matrix(y$x, nrow = p, ncol = p)
  y
}
