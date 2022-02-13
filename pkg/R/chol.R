## ID: chol.R, last updated 2022-02-12, F.Osorio

cholupdate <- function(r, x)
{ ## rank-1 update to Cholesky factorization
  if (is.data.frame(r))
    r <- as.matrix(r)
  if (!is.matrix(r))
    stop("supply a matrix-like 'r'")
  if (!is.numeric(r))
    stop("argument r is not a numeric matrix")

  dr <- dim(r)
  n <- dr[1]
  p <- dr[2]
  if (n != p)
    stop("argument r is not a square matrix")
  if (!is.upper.tri(r, diag = TRUE))
    stop("r must be an upper triangular matrix")
  storage.mode(r) <- "double"

  if (!is.vector(x))
    x <- as.vector(x)
  if (!is.numeric(x))
    stop("argument x is not a numeric vector")
  if (p != length(x))
    stop("arguments r and x are not compatible")

  z <- .C("chol_update",
          r = r,
          ldr = as.integer(n),
          p = as.integer(p),
          x = as.double(x))$r
  # output
  z[row(z) > col(z)] <- 0
  z
}
