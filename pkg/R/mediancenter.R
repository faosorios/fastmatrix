## ID: mediancenter.R, last updated 2021-10-16, F.Osorio

mediancenter <- function(x)
{ ## computes the mediancenter for a sample of multivariate observations
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  storage.mode(x) <- "double"
  z <- .Fortran("median_center",
                x = x,
                ldx = as.integer(n),
                n = as.integer(n),
                p = as.integer(p),
                median = double(p),
                iter = as.integer(0),
                info = as.integer(0))[c("median","iter")]
  z
}
