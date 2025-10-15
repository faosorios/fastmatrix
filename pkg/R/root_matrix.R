## ID: root_matrix.R, last updated 2025-10-14, F.Osorio

matrix.sqrt <- function(a, method = "DB", maxiter = 50, tol = 1e-8)
{ ## square root of a square matrix
  if (is.data.frame(a))
    a <- as.matrix(a)
  if (!is.matrix(a))
    stop("supply a matrix-like 'a'")
  if (!is.numeric(a))
    stop("argument a is not a numeric matrix")

  da <- dim(a)
  n <- da[1]
  p <- da[2]
  if (n != p)
    stop("argument a is not a square matrix")
  storage.mode(a) <- "double"
  method <- pmatch(method, c("DB", "schur"))

   switch(method,
         "DB" = {
            z <- .C("sqrt_mat_DB",
                    a = a,
                    lda = as.integer(n),
                    n = as.integer(n),
                    maxiter = as.integer(maxiter),
                    tol = as.double(tol),
                    iterations = as.integer(0))[c("a","iterations")]
            iterations <- z$iterations
            z <- z$a
            attr(z, 'iterations') <- iterations
         },
         "schur" = {
            z <- .C("sqrt_mat_schur",
                    a = a,
                    lda = as.integer(n),
                    n = as.integer(n))
            z <- z$a
         })
  z
}
