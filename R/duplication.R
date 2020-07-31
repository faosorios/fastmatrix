## ID: duplication.R, last updated 2020-07-31, F.Osorio

dupl.info <- function(n = 1)
{ ## returns compact information to form the duplication matrix of order 'n'
  n <- as.integer(n)
  row <- seq(n^2)
  col <- integer(n^2)
  z <- .C("dupl_cols",
    order = as.integer(n),
    col = as.integer(col))
  z <- list(row = row, col = z$col, order = n)
  z
}

duplication <- function(n = 1, matrix = FALSE)
{ ## returns the duplication matrix of order 'n'
  ## based on the post by Charles Berry, 2006-09-09 at R-help
  do.matrix <- matrix
  z <- dupl.info(n)
  if (do.matrix) {
    mat <- matrix(0, nrow = n^2, ncol = n * (n + 1) / 2)
    storage.mode(mat) <- "integer"
    col <- z$col
    z <- .Fortran("dupl_mat",
      mat = mat,
      ldmat = as.integer(n^2),
      order = as.integer(n),
      col = as.integer(col))$mat
  }
  z
}

dupl.prod <- function(n = 1, x)
{ ## returns y <- duplication(n) %*% x
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.vector(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  storage.mode(x) <- "double"

  rows <- n^2
  cols <- n * (n + 1) / 2
  dx <- dim(x)
  xrow <- dx[1]
  xcol <- dx[2]
  if (xrow != cols)
    stop("incompatible dimensions.")

  y <- matrix(0, nrow = rows, ncol = xcol)
  storage.mode(y) <- "double"
  col <- dupl.info(n)$col

  z <- .Fortran("dlmm",
    col = as.integer(col),
    order = as.integer(n),
    x = x,
    ldx = as.integer(xrow),
    xrow = as.integer(xrow),
    xcol = as.integer(xcol),
    y = y,
    ldy = as.integer(rows))
  z$y
}
