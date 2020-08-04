## ID: duplication.R, last updated 2020-08-01, F.Osorio

dupl.info <- function(n = 1, condensed = FALSE)
{ ## returns compact information to form the duplication matrix of order 'n'
  n <- as.integer(n)
  col <- integer(n^2)
  z <- .C("dupl_cols",
    order = as.integer(n),
    col = as.integer(col))
  z <- list(col = z$col, order = n)
  if (!condensed) {
    z$row <- seq(n^2)
    z <- z[c("row","col","order")]
  }
  z
}

duplication <- function(n = 1, matrix = FALSE, condensed = TRUE)
{ ## returns the duplication matrix of order 'n'
  ## based on the post by Charles Berry, 2006-09-09 at R-help
  do.matrix <- matrix
  z <- dupl.info(n, condensed = condensed)
  if (do.matrix) {
    mat <- matrix(0, nrow = n^2, ncol = n * (n + 1) / 2)
    storage.mode(mat) <- "integer"
    col <- z$col
    z <- .C("duplication_mat",
      mat = mat,
      ldmat = as.integer(n^2),
      order = as.integer(n),
      col = as.integer(col))$mat
  }
  z
}

dupl.prod <- function(n = 1, x, transposed = FALSE, side = "left")
{ ## returns y <- duplication(n) %*% x
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.vector(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  storage.mode(x) <- "double"

  dx <- dim(x)
  xrow <- dx[1]
  xcol <- dx[2]

  switch(side,
         "left" = {
           if (transposed) {
             rows <- n * (n + 1) / 2
             cols <- n^2
             if (xrow != cols)
               stop("incompatible dimensions.")
             y <- matrix(0, nrow = rows, ncol = xcol)
             storage.mode(y) <- "double"
             col <- dupl.info(n)$col
             counts <- table(col)
             oc <- order(col)
             col <- seq(n^2)[oc]
             z <- .C("dupl_left_trans",
                     x = x,
                     ldx  = as.integer(xrow),
                     xrow = as.integer(xrow),
                     xcol = as.integer(xcol),
                     col  = as.integer(col),
                     order = as.integer(n),
                     counts = as.integer(counts),
                     y = y,
                     ldy  = as.integer(rows))$y
              z
           } else {
             rows <- n^2
             cols <- n * (n + 1) / 2
             if (xrow != cols)
               stop("incompatible dimensions.")
             y <- matrix(0, nrow = rows, ncol = xcol)
             storage.mode(y) <- "double"
             col <- dupl.info(n)$col
             z <- .C("dupl_left_mult",
                     x = x,
                     ldx  = as.integer(xrow),
                     xrow = as.integer(xrow),
                     xcol = as.integer(xcol),
                     col  = as.integer(col),
                     order = as.integer(n),
                     y = y,
                     ldy  = as.integer(rows))$y
              z
           }
         },
         "right" = {
           if (transposed) {
             z <- NULL
           } else {
             z <- NULL
           }
         })
  z
}
