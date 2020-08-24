## ID: equilibrate.R, last updated 2020-08-22, F.Osorio

equilibrate <- function(x, scale = TRUE)
{ ## columns equilibration of a rectangular matrix to reduce its condition number
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix" )

  dx <- dim(x)
  xrow <- dx[1]
  xcol <- dx[2]
  storage.mode(x) <- "double"

  z <- .C("equilibrate",
          x = x,
          ldx  = as.integer(xrow),
          nrow = as.integer(xrow),
          ncol = as.integer(xcol),
          scales = double(xcol),
          condition = as.double(0),
          job = as.integer(scale))[c("x","scales","condition")]
  attr(z$x, 'scales') <- z$scales
  attr(z$x, 'condition') <- z$condition
  z$x
}
