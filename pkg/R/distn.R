## ID: distn.R, last updated 2025-12-19, F.Osorio

# pdf, cdf, quantile and RNG functions for the chi distribution

dchi <- function(x, df = 1, log = FALSE)
{ # density of the chi distribution
  if (df <= 0.0)
    stop("degrees of freedom must be non-negative.")

  n <- length(x)
  ndf <- length(df)

  y <- .C("d_chi",
          n = as.integer(n),
          y = double(n),
          x = as.double(x),
          df = as.double(df),
          ndf = as.integer(ndf),
          give.log = as.integer(log))$y
  y
}

pchi <- function(q, df = 1, lower.tail = TRUE, log.p = FALSE)
{ # distribution function of the chi distribution
  if (df <= 0.0)
    stop("degrees of freedom must be non-negative.")

  n <- length(q)
  ndf <- length(df)

  y <- .C("p_chi",
          n = as.integer(n),
          y = double(n),
          q = as.double(q),
          df = as.double(df),
          ndf = as.integer(ndf),
          lower.tail = as.integer(lower.tail),
          log.p = as.integer(log.p))$y
  y
}

qchi <- function(p, df = 1, lower.tail = TRUE, log.p = FALSE)
{ # quantile function of the chi distribution
  if (df <= 0.0)
    stop("degrees of freedom must be non-negative.")

  n <- length(p)
  ndf <- length(df)

  y <- .C("q_chi",
          n = as.integer(n),
          y = double(n),
          p = as.double(p),
          df = as.double(df),
          ndf = as.integer(ndf),
          lower.tail = as.integer(lower.tail),
          log.p = as.integer(log.p))$y
  y
}

rchi <- function(n, df = 1)
{ # generation of chi random variables
  if (df <= 0.0)
    stop("degrees of freedom must be non-negative.")

  ndf <- length(df)
  # calling C code
  x <- .C("r_chi",
          n = as.integer(n),
          x = double(n),
          df = as.double(df),
          ndf = as.integer(ndf))$x
  x
}
