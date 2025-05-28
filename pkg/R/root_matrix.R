## ID: root_matrix.R, last updated 2025-05-26, F.Osorio

matrix.sqrt <- function(a, maxiter = 50, tol = 1e-8)
{ ## rank-1 update to Cholesky factorization
  if (is.data.frame(a))
    a <- as.matrix(a)
  if (!is.matrix(a))
    stop("supply a matrix-like 'a'")
  if (!is.numeric(a))
    stop("argument a is not a numeric matrix")

  da <- dim(a)
  n <- da[1]
  p <- da[2]
  if (n != p)
    stop("argument a is not a square matrix")
  storage.mode(a) <- "double"

  z <- .C("sqrt_mat_DB",
          a = a,
          lda = as.integer(n),
          n = as.integer(n),
          info = as.integer(0),
          maxiter = as.integer(maxiter),
          tol = as.double(tol),
          iterations = as.integer(0))[c("a","iterations")]
  iterations <- z$iterations
  z <- z$a
  attr(z, 'iterations') <- iterations
  z
}
