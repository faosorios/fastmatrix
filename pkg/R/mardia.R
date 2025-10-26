## ID: mardia.R, last updated 2020-11-22, F.Osorio

kurtosis <- function(x)
{ ## Mardia's multivariate kurtosis coefficient
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  z <- cov.weighted(x, weights = rep(1, n))
  storage.mode(x) <- "double"
  storage.mode(z$cov) <- "double"

  z <- .C("skewness_and_kurtosis",
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          center = as.double(z$mean),
          cov = z$cov,
          stats = double(2),
          task = as.integer(0))$stats[2]
  attr(z, 'excess') <- z - p * (p + 2)
  z
}

skewness <- function(x)
{ ## Mardia's multivariate skewness coefficient
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  z <- cov.weighted(x, weights = rep(1, n))
  storage.mode(x) <- "double"
  storage.mode(z$cov) <- "double"

  z <- .C("skewness_and_kurtosis",
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          center = as.double(z$mean),
          cov = z$cov,
          stats = double(2),
          task = as.integer(1))$stats[1]
  z
}

mardia.test <-
function(x)
{ ## Mardia test for multivariate normality
  ## Mardia, K.V. (1970). Biometrika 57, 519-530.
  ## doi: 10.1093/biomet/57.3.519
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  z <- cov.weighted(x, weights = rep(1, n))
  storage.mode(x) <- "double"
  storage.mode(z$cov) <- "double"

  z <- .C("mardia_stat",
          x = x,
          n = as.integer(n),
          p = as.integer(p),
          center = as.double(z$mean),
          cov = z$cov,
          coef = double(2),
          stats = double(2))[c("coef","stats")]

  skewness <- z$coef[1]
  kurtosis <- z$coef[2]
  a  <- z$stats[1]
  b  <- z$stats[2]
  names(a) <- "Skewness"
  names(b) <- "Kurtosis"
  df <- p * (p + 1) * (p + 2) / 6
  a.pval <- 1 - pchisq(a, df = df)
  b.pval <- 2 * (1 - pnorm(abs(b)))
  method <- "Mardia test"

  ## output object
  z <- list(skewness = list(statistic = a, parameter = df, p.value = a.pval, coef = skewness),
            kurtosis = list(statistic = b, p.value = b.pval, coef = kurtosis), method = method)
  class(z) <- "Mardia.test"
  z
}

print.Mardia.test <- function(x, digits = 4, ...)
{
  cat("\n")
  cat(paste(x$method, "for multivariate normality", sep = " "), "\n")
  cat("\n")
  # 1st row: 
  row <- character()
  row <- c(row, paste(names(x$skewness$statistic), "statistic =",
                      format(round(x$skewness$statistic, digits = digits))))
  row <- c(row, paste("df =", x$skewness$parameter))
  row <- c(row, paste("p-value =", format(round(x$skewness$p.value, digits = digits))))
  cat(strwrap(paste(row, collapse = ", ")), sep = "\n")
  # 2nd row: 
  row <- character()
  row <- c(row, paste(names(x$kurtosis$statistic), "statistic =",
                      format(round(x$kurtosis$statistic, digits = digits))))
  row <- c(row, paste("p-value =", format(round(x$kurtosis$p.value, digits = digits))))
  cat(strwrap(paste(row, collapse = ", ")), sep = "\n")
  cat("alternative hypothesis: true distribution is not multivariate gaussian.\n")
  cat("\n")
  cat("sample skewness:", format(round(x$skewness$coef, digits = digits)))
  cat("\n")
  cat("sample kurtosis:", format(round(x$kurtosis$coef, digits = digits)))
  cat("\n")
  invisible(x)
}

