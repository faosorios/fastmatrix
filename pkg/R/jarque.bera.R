## ID: jarque.bera.R, last updated 2023-07-23, F.Osorio

JarqueBera.test <-
function(x, test = "DH")
{
  if (!is.vector(x))
    stop("argument 'x' must be a vector.")
  if (!is.numeric(x))
    stop("argument 'x' is not a numeric vector.")

  # removing NA's
  ok <- !is.na(x)
  x <- x[ok]

  n <- length(x)
  switch(test,
    "DH" = {
      if (n <= 8)
        stop("sample size must be greater than 8.")
      z <- .C("doornik_hansen", 
              x = as.double(x),
              n = as.integer(n),
              skewness = double(1),
              kurtosis = double(1),
              stat = double(1))[c("skewness","kurtosis","stat")]
      skewness <- z$skewness
      kurtosis <- z$kurtosis
      stat <- z$stat
      names(stat) <- "Doornik-Hansen"
      method <- "Doornik-Hansen omnibus test"
    },
    "JB" = {
      z <- .C("jarque_bera", 
              x = as.double(x),
              n = as.integer(n),
              skewness = double(1),
              kurtosis = double(1),
              stat = double(1))[c("skewness","kurtosis","stat")]
      skewness <- z$skewness
      kurtosis <- z$kurtosis
      stat <- z$stat
      names(stat) <- "Jarque-Bera"
      method <- "Jarque-Bera test"
    },
    "robust" = {
      u <- x - median(x)
      z <- .C("robust_JB", 
              x = as.double(x),
              u = as.double(u),
              n = as.integer(n),
              skewness = double(1),
              kurtosis = double(1),
              stat = double(1))[c("skewness","kurtosis","stat")]
      skewness <- z$skewness
      kurtosis <- z$kurtosis
      stat <- z$stat
      names(stat) <- "Robust Jarque-Bera"
      method <- "Robust Jarque-Bera test"
    },
    "ALM" = {
      z <- .C("urzua_ALM", 
              x = as.double(x),
              n = as.integer(n),
              skewness = double(1),
              kurtosis = double(1),
              stat = double(1))[c("skewness","kurtosis","stat")]
      skewness <- z$skewness
      kurtosis <- z$kurtosis
      stat <- z$stat
      names(stat) <- "adjusted LM"
      method <- "Adjusted Lagrange multiplier test"
    },
    stop(paste("unimplemented option:", test))
  )
  df <- 2
  pval <- 1 - pchisq(stat, df = df)

  ## output object
  z <- list(statistic = stat, parameter = df, p.value = pval, skewness = skewness, kurtosis = kurtosis, method = method)
  class(z) <- "JarqueBera.test"
  z
}

print.JarqueBera.test <- function(x, digits = 4, ...)
{
  cat("\n")
  cat(paste(x$method, "for normality", sep = " "), "\n")
  cat("\n")
  out <- character()
  out <- c(out, paste(names(x$statistic), "statistic =",
                      format(round(x$statistic, digits = digits))))
  out <- c(out, paste("df =", x$parameter))
  out <- c(out, paste("p-value =", format(round(x$p.value, digits = digits))))
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("alternative hypothesis: true distribution is not gaussian.\n")
  cat("\n")
  cat("sample skewness:", format(round(x$skewness, digits = digits)))
  cat("\n")
  cat("sample kurtosis:", format(round(x$kurtosis, digits = digits)))
  cat("\n")
  invisible(x)
}
