## ID: norms.R, last updated 2020-08-16, F.Osorio

matrix.norm <- function(x, type = "Frobenius")
{ ## Computes a matrix norm
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'.")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix.")

  job <- switch(type,
                "Inf"       = 0,
                "1"         = 1,
                "Frobenius" = 2,
                "maximum"   = 3,
                stop("type not implemented."))

  dx <- dim(x)
  n <- dx[1]
  p <- dx[2]
  storage.mode(x) <- "double"

  z <- .C("matrix_norm",
          x = x,
          ldx  = as.integer(n),
          nrow = as.integer(n),
          ncol = as.integer(p),
          job  = as.integer(job),
          val  = as.double(0))$val
  z
}

matrix.inner <- function(x, y = x)
{ ## Computes the inner product between x and y
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'.")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix.")
  if (is.vector(x))
    x <- as.matrix(x)

  if (is.data.frame(y))
    y <- as.matrix(y)
  if (!is.matrix(y))
    stop("supply a matrix-like 'x'.")
  if (!is.numeric(y))
    stop("argument y is not a numeric matrix.")
  if (is.vector(y))
    y <- as.matrix(y)

  dx <- dim(x)
  dy <- dim(y)
  if (!all(dx == dy))
    stop( "arguments x and y do not have the same order.")

  storage.mode(x) <- "double"
  storage.mode(y) <- "double"

  z <- .Fortran("inner_frobenius",
                x = x,
                ldx  = as.integer(dx[1]),
                y = y,
                ldy  = as.integer(dy[1]),
                nrow = as.integer(dy[1]),
                ncol = as.integer(dy[2]),
                val  = as.double(0))$val
  z
}

minkowski <- function(x, p = 2)
{ ## computes a p-norm
  if (!is.vector(x))
    stop("argument x must be a vector.")
  if (!is.numeric(x))
    stop("argument x is not a numeric vector.")
  if (p < 1)
    stop("p must be greather or equal than 1.")

  n <- length(x)
  if (p == 1) {
    z <- .C("norm_one",
            x   = as.double(x),
            inc = as.integer(1),
            n   = as.integer(n),
            val = as.double(0))$val
  } else if (p == 2) {
    z <- .C("norm_two",
            x   = as.double(x),
            inc = as.integer(1),
            n   = as.integer(n),
            val = as.double(0))$val
  } else if (p == Inf) {
    z <- .C("norm_inf",
            x   = as.double(x),
            inc = as.integer(1),
            n   = as.integer(n),
            val = as.double(0))$val
  } else {
    z <- .C("norm_minkowski",
            x   = as.double(x),
            inc = as.integer(1),
            n   = as.integer(n),
            p   = as.double(p),
            val = as.double(0))$val
  }
  z
}
