## ID: bezier.R, last updated 2021-11-14, F.Osorio

bezier <- function(x, y, ngrid = 200)
{ ## smoothing via Bezier curve

  # validating arguments
  if (!is.vector(x))
    x <- as.vector(x)
  if (!is.numeric(x))
    stop("argument x is not a numeric vector")
  if (!is.vector(y))
    y <- as.vector(y)
  if (!is.numeric(y))
    stop("argument y is not a numeric vector")

  # remove all NAs
  OK <- complete.cases(x, y)
  x <- x[OK]
  y <- y[OK]

  n <- length(x)
  if (n != length(y))
    stop("arguments x and y do not have the same length.")

  grid <- seq(0, 1, length = ngrid)

  z <- .C("bezier_smoother",
          x = as.double(x),
          y = as.double(y),
          n = as.integer(n),
          grid = as.double(grid),
          ngrid = as.integer(ngrid),
          xgrid = double(ngrid),
          ygrid = double(ngrid))[c("xgrid","ygrid")]
  z
}
