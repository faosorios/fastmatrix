## ID: jacobi.R, last updated 2021-02-19, F.Osorio

jacobi <- function(a, b, start, maxiter = 200, tol = 1e-7)
{ ## Jacobi iterative method for solving linear equations
  if (is.data.frame(a))
    a <- as.matrix(a)
  if (!is.matrix(a))
    stop("supply a matrix-like 'a'")
  if (!is.numeric(a))
    stop("argument 'a' is not a numeric matrix" )

  da <- dim(a)
  n <- da[1]
  p <- da[2]
  if (n != p)
    stop("argument 'a' is not a square matrix")
  storage.mode(a) <- "double"

  if (is.matrix(b)) {
    b <- b[,1]
    warning("'b' is a matrix, using its first column instead")
  }
  b <- as.vector(b)
  storage.mode(b) <- "double"

  if (p != length(b))
    stop("'a' and 'b' are not compatible")

  if (missing(start)) {
    warning("'start' is not provided, using zero instead")
    start <- double(p)
  }

  x <- start
  if (is.matrix(x)) {
    x <- x[,1]
    warning("'start' is a matrix, using its first column instead")
  }
  x <- as.vector(x)
  storage.mode(x) <- "double"

  if (p != length(x))
    stop("'b' and 'start' must have the same length")

  z <- .C("jacobi_solver",
          a = a,
          lda = as.integer(n),
          p = as.integer(p),
          b = b,
          x = x,
          maxiter = as.integer(maxiter),
          tol = as.double(tol),
          iter = as.integer(0),
          info = as.integer(0))[c("x","iter","info")]
  if (z$info < 0) {
    stop(paste("argument ", -z$info, " had an illegal value.", sep = ""))
  } else if (z$info > 0) {
    stop("matrix 'a' has an zero diagonal element.")
  }
  if (z$iter >= maxiter)
    stop("maximum number of iterations exceeded.")

  x <- z$x
  attr(x, 'iterations') <- z$iter
  x
}
