## ID: rank1_update.R, last updated 2026-01-12, F.Osorio

rank1.update <- function(a, alpha, u, diag = FALSE)
{ ## Rank-1 undate: mat <- a + alpha * u %*% t(u)
  ## if diag = TRUE, mat <- diag(a) + alpha * u %*% t(u)
  if (is.data.frame(a))
    a <- as.matrix(a)
  if (!is.numeric(a))
    stop("'a' is not a numeric argument")
  if (diag) {
    job <- 0
    if (is.matrix(a))
      a <- diag(a)
    if (!is.vector(a))
      stop("supply a vector-like 'a'")
    n <- length(a)
  } else {
    job <- 1
    if (!is.matrix(a))
      stop("supply a matrix-like 'a'")
    da <- dim(a)
    n <- da[1]
    p <- da[2]
    if (n != p)
      stop("argument 'a' is not a square matrix")
  }

  storage.mode(a) <- "double"
  z <- .C("rank1_update",
          mat = double(n * n),
          ldmat = as.integer(n),
          n = as.integer(n),
          job = as.integer(job),
          a = a,
          alpha = as.double(alpha),
          u = as.double(u))$mat
  z <- matrix(z, nrow = n, ncol = n)
  z
}
