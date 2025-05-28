## ID: hankel.R, last updated 2025-05-26, F.Osorio

hankel <- function(x, y = NULL)
{ ## function returns the n-by-n Hankel matrix based on 'x' and 'y'
  if (!is.vector(x))
    stop("supply a vector-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric vector" )
  n <- length(x)

  if (is.null(y)) {
    y <- rep(0, n)
  } else {
    if (!is.vector(y))
      stop("supply a vector-like 'x'")
    if (!is.numeric(x))
      stop("argument x is not a numeric vector" )
    if (n != length(y))
      stop("arguments x and y do not have the same length.")
  }

  z <- .Fortran("hankel_mat",
                x = as.double(x),
                y = as.double(y),
                n = as.integer(n),
                mat = double(n * n),
                ldmat = as.integer(n),
                info = as.integer(0))
  if (z$info)
    stop(paste("hankel_mat gave error code", z$info))
  z <- matrix(z$mat, nrow = n, ncol = n)
  z
}
