## ID: bits.R, last updated 2020-08-13, F.Osorio

revbits <- function(x)
{ ## linear bit reversal
  if (length(x)) {
    x <- .Fortran("bitreversal",
                  x = as.integer(x),
                  m = as.integer(length(x)))$x
  }
  x
}
