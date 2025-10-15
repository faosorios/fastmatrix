## ID: RNG.R, last updated 2023-07-23, F.Osorio

rmnorm <-
function(n = 1, mean = rep(0, nrow(Sigma)), Sigma = diag(length(mean)))
{ ## multivariate normal random generation
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(mean) != nrow(Sigma))
    stop("mean and sigma have non-conforming size")
  p <- nrow(Sigma)

  y <- matrix(0, nrow = n, ncol = p)
  # call C code
  y <- .C("rng_mnorm",
          y = as.double(y),
          n = as.integer(n),
          p = as.integer(p),
          mean = as.double(mean),
          Sigma = as.double(Sigma))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}

rsphere <- 
function(n = 1, p = 2)
{ ## random vector generation uniformly on the unitary sphere
  if (n <= 0)
    stop("n must be a positive integer")
  if (p <= 1)
    stop("dimension must be > 1")

  y <- matrix(0, nrow = n, ncol = p)
  # call C code
  y <- .C("rng_sphere", 
          y = as.double(y),
          n = as.integer(n),
          p = as.integer(p))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}

rball <- 
function(n = 1, p = 2)
{ ## random vector generation uniformly in the unitary ball
  if (n <= 0)
    stop("n must be a positive integer")
  if (p <= 1)
    stop("dimension must be > 1")

  y <- matrix(0, nrow = n, ncol = p)
  # call C code
  y <- .C("rng_ball", 
          y = as.double(y),
          n = as.integer(n),
          p = as.integer(p))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}
