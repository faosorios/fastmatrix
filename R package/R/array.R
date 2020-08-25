## ID: array.R, last updated 2020-08-15, F.Osorio

array.mult <- function(a, b, x)
{ ## returns the array multiplication: y <- a %*% x %*% b
  ## with x an array
  if (is.data.frame(a))
    a <- as.matrix(a)
  if (!is.matrix(a))
    stop("supply a matrix-like 'a'")
  if (!is.numeric(a))
    stop("argument a is not a numeric matrix")

  if (is.data.frame(b))
    b <- as.matrix(b)
  if (!is.matrix(b))
    stop("supply a matrix-like 'b'")
  if (!is.numeric(b))
    stop("argument b is not a numeric matrix" )

  if (!is.array(x))
    stop("supply an array-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric array" )

  da <- dim(a)
  arow <- da[1]
  acol <- da[2]

  db <- dim(b)
  brow <- db[1]
  bcol <- db[2]

  dx <- dim(x)
  if (length(dx) != 3)
    stop("argument x is not a 3D array")
  m <- dx[1]
  n <- dx[2]
  p <- dx[3]

  if (acol != n)
    stop("incompatible dimensions.")

  if (brow != p)
    stop("incompatible dimensions.")

  y <- array(0, dim = c(m, arow, bcol))

  storage.mode(a) <- "double"
  storage.mode(b) <- "double"
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"

  z <- .Fortran("arraymult",
                a = a,
                lda  = as.integer(arow),
                arow = as.integer(arow),
                acol = as.integer(acol),
                b = b,
                ldb  = as.integer(brow),
                brow = as.integer(brow),
                bcol = as.integer(bcol),
                x = x,
                ldx  = as.integer(m),
                n = as.integer(m),
                y = y,
                ldy  = as.integer(m),
                info = as.integer(0))
  if (z$info)
    stop(paste("arraymult gave error code", z$info))
  z$y
}

bracket.prod <- function(a, x)
{ ## returns the bracket product between matrix 'a' and array 'x'
  if (is.data.frame(a))
    a <- as.matrix(a)
  if (!is.matrix(a))
    stop("supply a matrix-like 'a'")
  if (!is.numeric(a))
    stop("argument a is not a numeric matrix" )

  if (!is.array(x))
    stop("supply an array-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric array" )

  da <- dim(a)
  arow <- da[1]
  acol <- da[2]

  dx <- dim(x)
  if (length(dx) != 3)
    stop("argument x is not a 3D array")
  m <- dx[1]
  n <- dx[2]
  p <- dx[3]

  if (acol != m)
    stop("incompatible dimensions.")

  y <- array(rep(0, arow * n * p), dim = c(arow, n, p))

  storage.mode(a) <- "double"
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"

  z <- .Fortran("bracketprod",
                a = a,
                lda  = as.integer(arow),
                arow = as.integer(arow),
                acol = as.integer(acol),
                x = x,
                ldx  = as.integer(m),
                n = as.integer(n),
                p = as.integer(p),
                y = y,
                ldy  = as.integer(arow),
                info = as.integer(0))
  if (z$info)
    stop(paste("bracketprod gave error code", z$info))
  z$y
}
