## ID: ldl.R, last updated 2021-02-17, F.Osorio

ldl <- function(x)
{ ## LDL decomposition of a square matrix
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix" )

  dx <- dim(x)
  n <- dx[1]
  p <- dx[2]
  if (n != p)
    stop("argument x is not a square matrix")
  storage.mode(x) <- "double"

  z <- .Fortran("ldl_dcmp",
                x = x,
                ldx = as.integer(n),
                p = as.integer(p),
                d = double(p),
                info = as.integer(0))[c("x","d")]
  # output
  lower <- z$x
  lower[row(lower) <= col(lower)] <- 0
  diag(lower) <- 1
  z <- list(lower = lower, d = z$d)
  z
}
