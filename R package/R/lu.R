## ID: lu.R, last updated 2020-09-15, F.Osorio

is.lu <- function(x) is.list(x) && inherits(x, "lu")

lu <- function(x) UseMethod("lu")

lu.default <- function(x)
{ ## LU factorization of a square matrix
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
  storage.mode(x) <- "double"

  z <- .C("lu_dcmp",
          lu = x,
          ldx = as.integer(n),
          n = as.integer(n),
          p = as.integer(p),
          pivot = integer(p))[c("lu","pivot")]
  class(z) <- "lu"
  z
}

extractL <- function(x)
{ # get L matrix from its LU factorization
  if (!is.lu(x)) stop("argument is not an LU factorization")
  L <- x$lu
  L[row(L) <= col(L)] <- 0
  diag(L) <- 1
  L
}

extractU <- function(x)
{ # get U matrix from its LU factorization
  if (!is.lu(x)) stop("argument is not an LU factorization")
  U <- x$lu
  U[row(U) > col(U)] <- 0
  U
}

constructX <- function(x)
{ # returns the original matrix from the 'LU' object
  if (!is.lu(x)) stop("argument is not an LU factorization")
  L <- extractL(x)
  ldx <- p <- nrow(L)
  storage.mode(L) <- "double"
  L <- .Fortran("pivot_mat",
                L = t(L),
                ldx = as.integer(ldx),
                p = as.integer(p),
                pivot = as.integer(x$pivot))$L
  U <- extractU(x)
  x <- crossprod(L, U)
  x
}

lu2inv <- function(x)
{ # computes the inverse of a matrix from its LU factorization
  if (!is.lu(x)) stop("argument is not an LU factorization")
  pivot <- x$pivot
  x <- x$lu
  storage.mode(x) <- "double"
  ldx <- p <- nrow(x)

  z <- .C("lu_inverse",
          x = x,
          ldx = as.integer(ldx),
          p = as.integer(p),
          pivot = as.integer(pivot))$x
  z
}

solve.lu <- function(a, b, ...)
{
  if(!is.lu(a))
    stop("this is the \"lu\" method for the generic function solve()")
  pivot <- a$pivot
  a <- a$lu
  storage.mode(a) <- "double"
  lda <- p <- nrow(a)

  if (missing(b)) {
    b <- diag(1, p)
    warning("b is missing use 'lu2inv' instead 'solve.lu' with an \"lu\" object")
  }
  b <- as.matrix(b)
  storage.mode(b) <- "double"
  nrhs <- ncol(b)

  if (p != nrow(b))
    stop("matrices a and b are not compatible")
  ldb <- p

  z <- .C("lu_solve",
          a = a,
          lda = as.integer(lda),
          p = as.integer(p),
          pivot = as.integer(pivot),
          b = b,
          ldb = as.integer(ldb),
          nrhs = as.integer(nrhs))$b
  if (nrhs == 1)
    z <- as.vector(z)
  z
}
