## ID: norms.R, last updated 2020-08-07, F.Osorio

matrix.norm <- function(x, type = "Frobenius")
{ ## Computes a matrix norm
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'.")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix.")

  job <- switch(type, "Frobenius" = "F",
                      "1"         = "1",
                      "inf"       = "I",
                      "maximum"   = "M")

  storage.mode(x) <- "double"
  z <- norm(x, type = job)
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
