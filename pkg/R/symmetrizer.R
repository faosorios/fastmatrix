## ID: symmetrizer.R, last updated 2020-09-16, F.Osorio

symm.info <- function(n = 1)
{ ## returns compact information to form the symmetrizer matrix of order 'n'
  n <- as.integer(n)
  row <- integer(n^2)
  # re-using commutation info
  row <- .Fortran("comm_rows",
                  m = as.integer(n),
                  n = as.integer(n),
                  row = as.integer(row))$row
  idx <- col <- seq(n^2)
  pos <- idx[row == col]

  row.nodiag <- row[-pos]
  col.nodiag <- col[-pos]
  val.nodiag <- rep(0.5, length(row.nodiag))

  val <- rep(0.5, n^2)
  row <- col <- idx
  val[pos] <- 1

  row <- c(row, row.nodiag)
  col <- c(col, col.nodiag)
  val <- c(val, val.nodiag)

  # output object
  z <- list(row = row, col = col, val = val, order = n)
  z
}

symmetrizer <- function(n = 1, matrix = FALSE)
{ ## wrapper to 'symm.info', if requested (but is not recommended) this
  ## function returns the symmetrizer matrix of order 'n', otherwise the
  ## compact information is returned.
  do.matrix <- matrix
  z <- symm.info(n)
  if (do.matrix) {
    mat <- matrix(0, nrow = n^2, ncol = n^2)
    storage.mode(mat) <- "double"
    len <- length(z$row)
    z <- .Fortran("symmetrizer_mat",
                  mat = mat,
                  ldmat = as.integer(n * n),
                  n = as.integer(n),
                  row = as.integer(z$row),
                  col = as.integer(z$col),
                  val = as.double(z$val),
                  len = as.integer(len),
                  info = as.integer(0))
    if (z$info)
      stop(paste("symmetrizer_mat gave error code", z$info))
    z <- z$mat
  }
  z
}

symm.prod <- function(n = 1, x = NULL, side = "left")
{ ## let N <- symmetrizer(n), be the symmetrizer matrix of order 'n'
  ## 'symm.prod' returns:
  ## y[,] <- N %*% x,    for side = "left",  or
  ## y[,] <- x %*% N,    for side = "right".
  if (is.null(x))
    stop("matrix 'x' is not supplied")
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.vector(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix")
  storage.mode(x) <- "double"

  dx <- dim(x)
  xrow <- dx[1]
  xcol <- dx[2]

  switch(side,
         "left" = { # y[,] <- N %*% x
           rows <- cols <- n^2
           if (cols != xrow)
             stop("incompatible dimensions.")
           y <- comm.prod(m = n, x = x, transposed = FALSE, side = "left")
           ldy <- nrow(y)
           storage.mode(y) <- "double"
           z <- .C("symmetrizer_prod",
                   x = 0.5 * x,
                   ldx = as.integer(xrow),
                   xrow = as.integer(xrow),
                   xcol = as.integer(xcol),
                   y = y,
                   ldy = as.integer(ldy))$x
            z
         },
         "right" = { # y[,] <- x %*% N
           rows <- cols <- n^2
           if (xcol != rows)
           stop("incompatible dimensions.")
           y <- comm.prod(m = n, x = x, transposed = FALSE, side = "right")
           ldy <- nrow(y)
           storage.mode(y) <- "double"
           z <- .C("symmetrizer_prod",
                   x = 0.5 * x,
                   ldx = as.integer(xrow),
                   xrow = as.integer(xrow),
                   xcol = as.integer(xcol),
                   y = y,
                   ldy = as.integer(ldy))$x
            z
         })
  z
}
