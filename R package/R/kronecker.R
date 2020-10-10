## ID: kronecker.R, last updated 2020-08-16, F.Osorio

kronecker.prod <- function(x, y = x)
{ ## returns the kronecker product of two matrices
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.vector(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix" )

  if (is.data.frame(y))
    y <- as.matrix(y)
  if (is.vector(y))
    y <- as.matrix(y)
  if (!is.matrix(y))
    stop("supply a matrix-like 'y'")
  if (!is.numeric(y))
    stop("argument y is not a numeric matrix" )

  dx <- dim(x)
  dy <- dim(y)
  z <- double(prod(dx * dy))
  dim(z) <- dx * dy

  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
  storage.mode(z) <- "double"

  z <- .C("kronecker_prod",
          x = x,
          xrow = as.integer(dx[1]),
          xcol = as.integer(dx[2]),
          y = y,
          yrow = as.integer(dy[1]),
          ycol = as.integer(dy[2]),
          z = z)$z
  z
}
