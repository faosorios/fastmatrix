## ID: asSymmetric.R, last updated 2021-04-07, F.Osorio

asSymmetric <- function(x, lower = TRUE)
{ ## coerces to symmetry
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument 'x' is not a numeric matrix")
  if (nrow(x) != ncol(x))
    stop("'x' must be a square matrix.")

  if (lower)
    x[row(x) < col(x)] <- x[row(x) > col(x)]
  else
    x[row(x) > col(x)] <- x[row(x) < col(x)]
  x
}
