## ID: utils.R, last updated 2020-08-10, F.Osorio

vec <- function(x)
{ ## returns a vector obtained by stacking the columns of x
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix" )

  if (is.vector(x))
    return(x) # quick return

  y <- as.double(x)
  y
}

vech <- function(x)
{ ## returns a vector obtained by stacking the lower triangular part of x
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix" )

  dx <- dim(x)
  xrow <- dx[1]
  xcol <- dx[2]
  if (xrow != xcol)
    stop("argument x is not a square matrix")

  n <- xcol * (xcol + 1) / 2
  y <- .C("mat2vech",
          x = as.double(x),
          ldx = as.integer(xrow),
          p = as.integer(xcol),
          y = double(n))$y
  y
}
