## ID: cg.R, last updated 2021-02-20, F.Osorio

cg <- function(a, b, maxiter = 200, tol = 1e-7)
{ ## conjugate gradients iterative method for solving linear equations
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

  z <- .C("cg_solver",
          a = a,
          lda = as.integer(n),
          p = as.integer(p),
          b = b,
          x = double(p),
          maxiter = as.integer(maxiter),
          tol = as.double(tol),
          iter = as.integer(0),
          info = as.integer(0))[c("x","iter","info")]
  if (z$info < 0)
    stop(paste("argument ", -z$info, " had an illegal value.", sep = ""))
  if (z$iter >= maxiter)
    stop("maximum number of iterations exceeded.")

  x <- z$x
  attr(x, 'iterations') <- z$iter
  x
}
