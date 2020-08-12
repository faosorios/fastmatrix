## ID: norms.R, last updated 2020-08-07, F.Osorio

matrix.norm <- function(x, type = "Frobenius")
{ ## Computes a matrix norm
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix" )

  storage.mode(x) <- "double"
  dx <- dim(x)

  switch(type,
         "1" = {
           z <- .Fortran("maxcol_norm",
                         x = x,
                         ldx = as.integer(dx[1]),
                         n = as.integer(dx[1]),
                         p = as.integer(dx[2]),
                         val = as.double(0))$val
          },
          "Frobenius" = {
            n <- prod(dx)
            z <- .C("norm_two",
                    x = x,
                    inc = as.integer(1),
                    n = as.integer(n),
                    val = as.double(0))$val
          },
          "inf" = {
            z <- .Fortran("maxrow_norm",
                          x = x,
                          ldx = as.integer(dx[1]),
                          n = as.integer(dx[1]),
                          p = as.integer(dx[2]),
                          val = as.double(0))$val
          })
  z
}
