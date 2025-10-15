## ID: condition.R, last updated 2023-10-10, F.Osorio

scaled.condition <- function(x, scales = FALSE)
{ ## scaled condition number
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument 'x' is not a numeric matrix")
  if (nrow(x) == ncol(x))
    stop("'x' must be a rectangular matrix.")

  z <- equilibrate(x, scale = TRUE)
  colScales <- 1. / attr(z, 'scales')
  d <- svd(z, nu = 0, nv = 0)$d
  p <- length(d)
  cn <- d[1] / d[p]
  if (scales)
    attr(cn, 'scales') <- colScales
  cn
}
