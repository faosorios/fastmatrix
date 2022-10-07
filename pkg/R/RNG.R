## ID: RNG.R, last updated 2022-08-19, F.Osorio

rmnorm <-
function(n = 1, mean = rep(0, nrow(Sigma)), Sigma = diag(length(mean)))
{ # multivariate normal random generation
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(mean) != nrow(Sigma))
    stop("mean and sigma have non-conforming size")
  p <- nrow(Sigma)

  y <- matrix(0, nrow = n, ncol = p)
  # call C code
  y <- .C("mnorm_rand",
          y = as.double(y),
          n = as.integer(n),
          p = as.integer(p),
          mean = as.double(mean),
          Sigma = as.double(Sigma))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}
