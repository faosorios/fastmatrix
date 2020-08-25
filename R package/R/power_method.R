## ID: power_method.R, last updated 2020-08-22, F.Osorio

power.method <- function(x, only.value = FALSE, maxiter = 100, tol = 1e-8)
{ ## power method to approximate dominant eigenvalue and eigenvector
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix" )

  dx <- dim(x)
  n <- dx[1]
  p <- dx[2]
  if (n != p)
    stop("argument x is not a square matrix")
  if (!isSymmetric(x))
    stop("only implemented for symmetric matrices.")

  storage.mode(x) <- "double"

  # initial estimate of 1st eigenvector
  vector <- double(p)
  vector[1] <- 1.0

  z <- .C("power_method",
          x = x,
          ldx  = as.integer(n),
          n = as.integer(n),
          p = as.integer(p),
          vector = as.double(vector),
          value = as.double(0),
          maxiter = as.integer(maxiter),
          tol = as.double(tol),
          numIter = as.integer(0))[c("value","vector","numIter")]
  numIter <- z$numIter
  if (only.value)
    z <- z$value
  else
    z <- z[c("value","vector")]
  attr(z, 'iterations') <- numIter
  z
}
