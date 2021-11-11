## ID: equilibrate.R, last updated 2021-11-10, F.Osorio

equilibrate <- function(x, scale = TRUE)
{ ## equilibration of a rectangular (only columns) or a symmetric matrix to reduce
  ## its condition number
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

  task <- isSymmetric(unname(x))

  if (task) {
    symbolname <- "equilibrate_cols"
    z <- .C("equilibrate_sym",
            x = x,
            ldx = as.integer(xrow),
            p = as.integer(xcol),
            scales = double(xcol),
            condition = as.double(0),
            largest = as.double(0),
            info = as.integer(0))[c("x","scales","condition","info")]
    if (scale) { # equilibrating if requested
      z$x <- .C("equilibrating_sym",
                x = x,
                p = as.integer(xcol),
                scales = as.double(z$scales))$x
    }
  } else {
    symbolname <- "equilibrate_cols"
    z <- .Fortran("equilibrate_cols",
                  x = x,
                  ldx = as.integer(xrow),
                  nrow = as.integer(xrow),
                  ncol = as.integer(xcol),
                  scales = double(xcol),
                  condition = as.double(0),
                  job = as.integer(scale),
                  info = as.integer(0))[c("x","scales","condition","info")]
  }

  if (z$info)
    stop(paste(symbolname, "gave error code", z$info))
  attr(z$x, 'scales') <- z$scales
  attr(z$x, 'condition') <- 1. / z$condition
  z$x
}
