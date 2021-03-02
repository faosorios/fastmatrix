## ID: commutation.R, last updated 2020-08-26, F.Osorio

comm.info <- function(m = 1, n = m, condensed = TRUE)
{ ## returns compact information to form the commutation matrix of order 'mn'
  m <- as.integer(m)
  n <- as.integer(n)
  row <- integer(m * n)
  z <- .Fortran("comm_rows",
                m = as.integer(m),
                n = as.integer(n),
                row = as.integer(row))
  z <- list(row = z$row, m = m, n = n)
  if (!condensed) {
    z$col <- seq(m * n)
    z <- z[c("row","col","m","n")]
  }
  z
}

commutation <- function(m = 1, n = m, matrix = FALSE, condensed = FALSE)
{ ## wrapper to 'comm.info', if requested (but is not recommended) this
  ## function returns the commutation matrix of order 'mn', otherwise the
  ## compact information is returned.
  ## based on the code posted at:
  ## https://math.stackexchange.com/questions/307299/kronecker-product-and-the-commutation-matrix
  do.matrix <- matrix
  z <- comm.info(m, n, condensed = condensed)
  if (do.matrix) {
    mat <- matrix(0, nrow = m * n, ncol = m * n)
    storage.mode(mat) <- "integer"
    row <- z$row
    z <- .Fortran("commutation_mat",
                  mat = mat,
                  ldmat = as.integer(m * n),
                  m = as.integer(m),
                  n = as.integer(n),
                  row = as.integer(row),
                  info = as.integer(0))
    if (z$info)
      stop(paste("commutation_mat gave error code", z$info))
    z <- z$mat
  }
  z
}

comm.prod <- function(m = 1, n = m, x = NULL, transposed = FALSE, side = "left")
{ ## let K <- commutation(m, n), be the commutation matrix of order 'mn'
  ## 'comm.prod' returns:
  ## y[,] <- K %*% x,    for transposed = FALSE, side = "left",  or
  ## y[,] <- t(K) %*% x, for transposed = TRUE,  side = "left",  or
  ## y[,] <- x %*% K,    for transposed = FALSE, side = "right", or
  ## y[,] <- x %*% t(K), for transposed = TRUE,  side = "right".
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
         "left" = {
           if (transposed) { # y[,] <- t(K) %*% x
             rows <- cols <- m * n
             if (cols != xrow)
               stop("incompatible dimensions.")
             y <- matrix(0, nrow = rows, ncol = xcol)
             storage.mode(y) <- "double"
             col <- comm.info(n, m)$row
             col <- order(col)
             z <- .Fortran("comm_left_mult",
                           col = as.integer(col),
                           m = as.integer(m),
                           n = as.integer(n),
                           x = x,
                           ldx  = as.integer(xrow),
                           xrow = as.integer(xrow),
                           xcol = as.integer(xcol),
                           y = y,
                           ldy  = as.integer(rows),
                           info = as.integer(0))[c("y","info")]
              if (z$info)
                stop(paste("comm_left_mult gave error code", z$info))
              z <- z$y
              z
           } else { # y[,] <- K %*% x
             rows <- cols <- m * n
             if (cols != xrow)
               stop("incompatible dimensions.")
             y <- matrix(0, nrow = rows, ncol = xcol)
             storage.mode(y) <- "double"
             col <- comm.info(m, n)$row
             col <- order(col)
             z <- .Fortran("comm_left_mult",
                           col = as.integer(col),
                           m = as.integer(m),
                           n = as.integer(n),
                           x = x,
                           ldx  = as.integer(xrow),
                           xrow = as.integer(xrow),
                           xcol = as.integer(xcol),
                           y = y,
                           ldy  = as.integer(rows),
                           info = as.integer(0))[c("y","info")]
              if (z$info)
                stop(paste("comm_left_mult gave error code", z$info))
              z <- z$y
              z
           }
         },
         "right" = {
           if (transposed) { # y[,] <- x %*% t(K)
             rows <- cols <- m * n
             if (xcol != rows)
               stop("incompatible dimensions.")
             y <- matrix(0, nrow = xrow, ncol = cols)
             storage.mode(y) <- "double"
             col <- comm.info(n, m)$row
             z <- .Fortran("comm_right_mult",
                           col = as.integer(col),
                           m = as.integer(m),
                           n = as.integer(n),
                           x = x,
                           ldx  = as.integer(xrow),
                           xrow = as.integer(xrow),
                           xcol = as.integer(xcol),
                           y = y,
                           ldy  = as.integer(xrow),
                           info = as.integer(0))[c("y","info")]
              if (z$info)
                stop(paste("comm_right_mult gave error code", z$info))
              z <- z$y
              z
           } else { # y[,] <- x %*% K
             rows <- cols <- m * n
             if (xcol != rows)
               stop("incompatible dimensions.")
             y <- matrix(0, nrow = xrow, ncol = cols)
             storage.mode(y) <- "double"
             col <- comm.info(m, n)$row
             z <- .Fortran("comm_right_mult",
                           col = as.integer(col),
                           m = as.integer(m),
                           n = as.integer(n),
                           x = x,
                           ldx  = as.integer(xrow),
                           xrow = as.integer(xrow),
                           xcol = as.integer(xcol),
                           y = y,
                           ldy  = as.integer(xrow),
                           info = as.integer(0))[c("y","info")]
              if (z$info)
                stop(paste("comm_right_mult gave error code", z$info))
              z <- z$y
              z
           }
         })
  z
}
