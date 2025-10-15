## ID: fnc_matrix.R, last updated 2025-10-12, F.Osorio

matrix.fun <- function(a, FUN = "log")
{ ## Parlett method for a function of a triangular matrix
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
  if (!is.upper.tri(a))
    stop("argument a is not an upper triangular matrix")
  storage.mode(a) <- "double"

  # evaluating diagonal elements
  f <- matrix(0, nrow = n, ncol = n)
  b <- diag(a)
  FUN <- match.fun(FUN)
  b <- FUN(b)
  bad <- any(is.infinite(b)) || any(is.nan(b)) || anyNA(b)
  if (bad)
    stop("some illegal values were obtained")

  diag(f) <- b
  storage.mode(f) <- "double"

  z <- .Fortran("fnc_parlett",
          a = a,
          lda = as.integer(n),
          n = as.integer(n),
          f = f,
          ldf = as.integer(n))
  z <- matrix(z$f, nrow = n, ncol = n)
  z
}
