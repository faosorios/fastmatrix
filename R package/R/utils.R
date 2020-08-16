## ID: utils.R, last updated 2020-08-10, F.Osorio

hadamard <- function(x, y)
{ ## returns the element-wise product of two matrices x and y
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix" )
  if (is.vector(x))
    x <- as.matrix(x)

  if (is.data.frame(y))
    y <- as.matrix(y)
  if (!is.matrix(y))
    stop("supply a matrix-like 'y'")
  if (!is.numeric(y))
    stop("argument y is not a numeric matrix" )
  if (is.vector(y))
    y <- as.matrix(y)

  dx <- dim(x)
  dy <- dim(y)
  if (!all(dx == dy))
    stop( "arguments x and y do not have the same order.")

  n <- prod(dx)
  z <- .Fortran("hadamard",
                x = as.double(x),
                y = as.double(y),
                n = as.integer(n),
                prod = double(n))$prod
  z <- matrix(z, nrow = dx[1], ncol = dx[2])
  z
}

sweep.operator <- function(x, k = 1, reverse = FALSE)
{ ## Gauss-Jordan sweep operator
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix" )
  if (is.vector(x))
    x <- as.matrix(x)

  dx <- dim(x)
  n <- dx[1]
  p <- dx[2]
  if (n != p)
    stop("argument x is not a square matrix")

  storage.mode(x) <- "double"
  z <- .Fortran("sweepop",
                x = x,
                ldx = as.integer(n),
                p = as.integer(p),
                k = as.integer(k),
                work = double(p),
                reverse = as.integer(reverse),
                info = as.integer(0))
  if (z$info)
    warning("error code", z$info, "in sweep.operator")

  z$x
}
