## ID: householder.R, last updated 2026-01-13, F.Osorio

house <- function(x, matrix = FALSE)
{ ## Householder transformation
  if (!is.vector(x))
    stop("supply a vector-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric vector")
  n <- length(x)

  z <- .C("house_vec",
          x = as.double(x),
          n = as.integer(n),
          u = double(n),
          tau = as.double(0))[c("u","tau")]
  if (matrix) {
    mat <- .C("house_mat",
              u = as.double(z$u),
              tau = as.double(z$tau),
              n = as.integer(n),
              mat = double(n * n))$mat
    mat <- matrix(mat, nrow = n, ncol = n)
  }

  o <- list(u = z$u, tau = z$tau)
  if (matrix)
    o <- mat
  o
}

house.prod <- function(a, x, side = "left")
{ ## apply a Householder reflection, 'house.prod' returns:
  ## y[,] <- H %*% a, side = "left", or
  ## y[,] <- a %*% H, side = "right",
  ## where H is the Householder transformation matrix.
  if (is.null(a))
    stop("argument 'x' is not supplied")
  if (is.data.frame(a))
    a <- as.matrix(a)
  if (!is.numeric(a))
    stop("'a' is not a numeric argument")
  if (is.vector(a)) {
    task <- 0
    n <- length(a)
  } else if (is.matrix(a)) {
    task <- 1
    da <- dim(a)
    n <- da[1]
    p <- da[2]
  } else 
    stop("the argument 'a' is not the right type.")

  if (!is.vector(x))
    stop("supply a vector-like 'x'")
  if (!is.numeric(x))
    stop("'x' is not a numeric argument")

  z <- house(x)
  if (task) {
    storage.mode(a) <- "double"
    switch(side,
           "left" = { # y[,] <- H %*% a
             if (n != length(x))
               stop("incompatible dimensions.")
             z <- .C("house_prod_mat",
                     a = a,
                     lda = as.integer(n),
                     nrow = as.integer(n),
                     ncol = as.integer(p),
                     job = as.integer(1),
                     u = as.double(z$u),
                     tau = as.double(z$tau))$a
             #z <- matrix(z, nrow = n, ncol = p)
           },
           "right" = { # y[,] <- a %*% H
             if (p != length(x))
               stop("incompatible dimensions.")
             z <- .C("house_prod_mat",
                     a = a,
                     lda = as.integer(n),
                     nrow = as.integer(n),
                     ncol = as.integer(p),
                     job = as.integer(0),
                     u = as.double(z$u),
                     tau = as.double(z$tau))$a 
           },
           stop("method not implemented."))
  } else {
    if (n != length(x))
      stop("incompatible dimensions.")
    z <- .C("house_prod_vec",
            x = as.double(x),
            n = as.integer(n),
            u = as.double(z$u),
            tau = as.double(z$tau))$x
  }
  z
}
