## ID: schur.R, last updated 2025-10-13, F.Osorio

schur <- function(x, vectors = TRUE) 
{ ## Schur decomposition of a nonsymmetric matrix
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix" )
  if (is.complex(x))
    stop("complex matrices not permitted at present")
  if (!all(is.finite(x))) 
    stop("infinite or missing values in 'x'")

  dx <- dim(x)
  n <- dx[1]
  p <- dx[2]
  if (n != p)
    stop("argument x is not a square matrix")
  storage.mode(x) <- "double"

  if (vectors) {
    v <- matrix(0, nrow = n, ncol = n)
    storage.mode(v) <- "double"
    ldv <- n
  } else {
    v <- 0 # dummy argument
    storage.mode(v) <- "double"
    ldv <- 1
  }

  z <- .C("schur_dcmp",
          x = x,
          ldx = as.integer(n),
          n = as.integer(n),
          task = as.integer(vectors),
          re = double(n),
          im = double(n),
          v = v,
          ldv = as.integer(ldv),
          info = as.integer(0))[c("x","re","im","v","info")]
  if (z$info < 0)
    stop(paste("argument ", -z$info, " had an illegal value.", sep = ""))
  if (z$info > 0)
    stop(paste("DGEES gave error code ", z$info, sep = ""))
  if (all(z$im == 0))
    values <- z$re
  else
    values <- complex(real = z$re, imaginary = z$im)
  # output
  o <- list(m = z$x, values = values)
  if (vectors)
    o$vectors <- z$v
  o
}
