## ID: krylov.R, last updated 2022-08-01, F.Osorio

krylov <- function(a, b, m = ncol(a))
{ ## construcs the Krylov matrix based on an n-by-n matrix a and an n-vector b

  ## validating arguments
  if (is.data.frame(a))
    a <- as.matrix(a)
  if (!is.matrix(a))
    stop("supply a matrix-like 'a'")
  if (!is.numeric(a))
    stop("argument 'a' is not a numeric matrix")
  da <- dim(a)
  n <- da[1]
  if (n != da[2])
    stop("argument 'a' is not a square matrix")
  storage.mode(a) <- "double"

  if (!is.vector(b))
    b <- as.vector(b)
  if (!is.numeric(b))
    stop("argument 'b' is not a numeric vector")
  if (n != length(b))
    stop("argument 'a' and 'b' are not compatible")

  if (m <= 0)
    stop("argument 'm' must be positive")

  mat <- matrix(0, nrow = n, ncol = m)
  storage.mode(mat) <- "double"
  z <- .C("krylov_mat",
          a = a,
          lda = as.integer(n),
          n = as.integer(n),
          b = as.double(b),
          m = as.integer(m),
          mat = mat,
          ldmat = as.integer(n),
          info = as.integer(0))
  if (z$info)
    stop(paste("krylov_mat gave error code", z$info))
  z <- z$mat
  z
}
