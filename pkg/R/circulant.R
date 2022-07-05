## ID: circulant.R, last updated 2022-02-26, F.Osorio

circulant <- function(x)
{ ## function returns the n-by-n circulant matrix based on 'x'
  if (!is.vector(x))
    stop("supply a vector-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric vector" )
  n <- length(x)
  z <- .Fortran("circulant_mat",
                x = as.double(x),
                n = as.integer(n),
                mat = double(n * n),
                ldmat = as.integer(n),
                info = as.integer(0))
  if (z$info)
    stop(paste("circulant_mat gave error code", z$info))
  z <- matrix(z$mat, nrow = n, ncol = n)
  z
}
