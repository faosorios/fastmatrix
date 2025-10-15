## ID: polymat.R, last updated 2022-08-08, F.Osorio

matrix.polynomial <- function(a, coef = rep(1, power + 1), power = length(coef))
{ ## evaluates a real general matrix polynomial using a Horner's scheme

  # validating arguments
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

  if (!is.vector(coef))
    coef <- as.vector(coef)
  if (!is.numeric(coef))
    stop("argument 'coef' is not a numeric vector")
  ncoef <- length(coef)

  if (all(coef == 1.0) && (power < 0))
    stop("argument 'power' must be positive")

  mat <- matrix(0, nrow = n, ncol = n)
  storage.mode(mat) <- "double"
  z <- .C("matrix_polynomial",
          a = a,
          lda = as.integer(n),
          n = as.integer(n),
          coef = as.double(coef),
          ncoef = as.integer(ncoef),
          mat = mat,
          ldmat = as.integer(n),
          info = as.integer(0))[c("mat", "info")]
  if (z$info)
    stop(paste("matrix_polynomial gave error code", z$info))
  z <- z$mat
  z
}
