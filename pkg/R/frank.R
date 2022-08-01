## ID: frank.R, last updated 2022-07-31, F.Osorio

frank <- function(n = 1)
{ ## function returns the n-by-n Frank matrix
  if (n <= 0)
    stop("'n' must be in non-negative")
  x <- matrix(0, nrow = n, ncol = n)
  storage.mode(x) <- "double"
  z <- .Fortran("frank_mat",
                x = x,
                ldx = as.integer(n),
                n = as.integer(n),
                info = as.integer(0))
  if (z$info)
    stop(paste("frank_mat gave error code", z$info))
  z <- z$x
  z
}
