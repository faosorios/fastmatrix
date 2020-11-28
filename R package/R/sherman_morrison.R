## ID: sherman_morrison.R, last updated 2020-11-18, F.Osorio

sherman.morrison <- function(a, b, d = b, inverted = FALSE)
{ ## Sherman-Morrison formula
  if (is.data.frame(a))
    a <- as.matrix(a)
  if (!is.matrix(a))
    stop("supply a matrix-like 'a'")
  if (!is.numeric(a))
    stop("argument a is not a numeric matrix" )

  da <- dim(a)
  n <- da[1]
  p <- da[2]
  if (n != p)
    stop("argument 'a' is not a square matrix")
  storage.mode(a) <- "double"

  z <- .C("sherman_morrison",
          a = a,
          lda  = as.integer(n),
          p = as.integer(p),
          b = as.double(b),
          d = as.double(d),
          inverted = as.integer(inverted))$a
  z
}
