## ID: wilson_hilferty.R, last updated 2020-11-24, F.Osorio

wilson.hilferty <- function(x)
{ # Wilson-Hilferty transformation
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  z <- cov.weighted(x, weights = rep(1, n))
  D2 <- Mahalanobis(x, z$mean, z$cov, inverted = FALSE)

  z <- .C("wilson_hilferty_chisq",
          distances = as.double(D2),
          n = as.integer(n),
          p = as.integer(p),
          z = double(n))$z
  z
}
