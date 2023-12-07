## ID: chol.R, last updated 2023-12-06, F.Osorio

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

mchol <- function(x)
{ ## modified Cholesky factorization of a symmetric but 
  ## not necessarily positive definite matrix
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix")

  dx <- dim(x)
  n <- dx[1]
  p <- dx[2]
  if (n != p)
    stop("argument x is not a square matrix")
  if (!isSymmetric(x))
    stop("argument x is not a symmetric matrix")
  storage.mode(x) <- "double"
  macheps <- .Machine$double.eps

  z <- .Fortran("mchol_dcmp",
                x = x,
                ldx = as.integer(n),
                p = as.integer(p),
                d = double(p),
                macheps = as.double(macheps),
                info = as.integer(0))[c("x","d","info")]  
  if (z$info != 0)
    warning(paste("mchol_dcmp gave error code", z$info))

  # output
  lower <- z$x
  lower[row(lower) <= col(lower)] <- 0
  diag(lower) <- 1
  z <- lower %*% diag(sqrt(z$d))
  z
}
