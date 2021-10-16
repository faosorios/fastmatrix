## ID: helmert.R, last updated 2021-05-09, F.Osorio

helmert <- function(n = 1)
{ ## function returns the Helmert matrix of order 'n'
  mat <- matrix(0, nrow = n, ncol = n)
  storage.mode(mat) <- "double"
  z <- .Fortran("helmert_mat",
                mat = mat,
                ldmat = as.integer(n),
                n = as.integer(n),
                info = as.integer(0))
  if (z$info)
    stop(paste("helmert_mat gave error code", z$info))
  z <- z$mat
  z
}
