## ID: sweep.R, last updated 2020-08-22, F.Osorio

sweep.operator <- function(x, k = 1, reverse = FALSE)
{ ## Gauss-Jordan sweep operator for symmetric matrices
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
  k <- as.vector(k)
  r <- length(k)

  storage.mode(x) <- "double"
  z <- .C("sweep_operator",
          x = x,
          ldx = as.integer(n),
          p = as.integer(p),
          k = as.integer(k),
          r = as.integer(r),
          reverse = as.integer(reverse))$x
  z
}
