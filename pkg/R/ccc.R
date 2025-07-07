## ID: ccc.R, last updated 2025-06-25, F.Osorio

ccc <- function(x, data, method = "z-transform", level = 0.95, equal.means = FALSE, ustat = TRUE, subset, na.action)
{ # estimation of the Lin's concordance correlation coefficient
  Call <- match.call()
  if (missing(x))
    stop("'x' is not supplied")
  if (inherits(x, "formula")) {
    mt <- terms(x, data = data)
    if (attr(mt, "response") > 0)
      stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    mf <- match.call(expand.dots = FALSE)
    names(mf)[names(mf) == "x"] <- "formula"
    mf$method <- mf$level <- mf$equal.means <- mf$ustat <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    na.act <- attr(mf, "na.action")
    z <- model.matrix(mt, mf)
  }
  else {
    z <- as.matrix(x)
    if (!missing(subset))
      z <- z[subset, , drop = FALSE]
    if (!missing(na.action))
      z <- na.omit(z)
    else
      z <- na.fail(z)
  }
  if (!is.numeric(z))
    stop("'ccc' applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]
  if (p > 2)
    stop("'ccc' is not implemented for p > 2")
  if ((level <= 0) || (level >= 1))
    stop("'level' must be in (0,1).")

  ## estimating mean vector and covariance matrix
  ones <- rep(1, n)
  storage.mode(z) <- "double"
  o <- .C("cov_weighted",
          x = z,
          n = as.integer(n),
          p = as.integer(p),
          weights = as.double(ones),
          mean = double(p),
          cov  = double(p * p))[c("mean","cov")]
  xcov <- matrix(o$cov, nrow = p, ncol = p)
  xbar <- o$mean
  
  ## restricted estimation
  if (equal.means) {
    lambda <- sum(solve(xcov, xbar)) / sum(solve(xcov))
    dev <- xbar - rep(lambda, p)
    Sigma <- xcov + outer(dev, dev)
  }

  ## computing CCC
  phi <- xcov[lower.tri(xcov, diag = TRUE)]
  diff  <- xbar[1] - xbar[2]
  ratio <- phi[1] / phi[3]
  a <- diff / ((phi[1] * phi[3])^.25) # location shift
  b <- sqrt(ratio) # scale shift
  rhoc <- 2. * phi[2] / (phi[1] + phi[3] + diff^2)
  accu <- 2. / (b + 1 / b + a^2)
  r <- rhoc / accu
  if (equal.means) {
    phi  <- Sigma[lower.tri(Sigma, diag = TRUE)]
    rho0 <- 2. * phi[2] / (phi[1] + phi[3])
    ratio <- phi[1] / phi[3]
    b0 <- sqrt(ratio) # scale shift
    acc0 <- 2. / (b0 + 1 / b0)
    r0 <- phi[2] / sqrt(phi[1] * phi[3])
    rel <- rho0 / r0
    var.rho0 <- (1 - r0^2) * (1 - rho0^2)
    var.rho0 <- var.rho0 * rel^2
    var.rho0 <- var.rho0 / (n - 2)
  }

  ## asymptotic variance using normal approximation (correction by Lin, Biometrics 56, pp.325, 2000)
  rel <- rhoc / r
  var.rhoc <- (1 - r^2) * (1 - rhoc^2) + 2 * rhoc * (1 - rhoc) * r * a^2 - 0.5 * rhoc^2 * a^4
  var.rhoc <- var.rhoc * rel^2
  var.rhoc <- var.rhoc / (n - 2)
  ## z-transformation
  var.z <- var.rhoc / (1 - rhoc^2)^2

  ## estimation using U-statistics (King & Chinchilli, J. Biopharm. Stat. 11, 83-105, 2001;
  ## Stat. Med. 20, 2131-2147, 2001)
  if (ustat)
    u <- ustat.rhoc(z)

  ## info for Bland & Altman plot
  average <- apply(z, 1, mean)
  diff <- z[,1] - z[,2]
  bland <- list(average = average, difference = diff)
  bland <- as.data.frame(bland)

  # confidence intervals for the concordance correlation coefficient
  z2rho <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)
  a  <- (1 - level) / 2
  a  <- c(a, 1 - a)
  qz <- qnorm(a)
  switch(method,
         "asymp" = { # asymptotic interval
           cname <- " CCC"
           cf <- rhoc
           SE <- sqrt(var.rhoc)
           ci <- cf + SE %o% qz
         },
         "z-transform" = { # z-transform (default)
           cname <- " z-trans"
           cf <- atanh(rhoc)
           SE <- sqrt(var.z)
           ci <- cf + SE %o% qz
           cf <- z2rho(cf)
           SE <- z2rho(SE)
           ci <- z2rho(ci)
         },
         stop("method = ", method, " is not implemented."))
  ci <- c(cf, SE, ci)
  names(ci) <- c(cname, "SE", "lower", "upper")

  ## creating the output object
  out <- list(call = Call, x = z, dims = dz, ccc = rhoc, var.ccc = var.rhoc, accuracy = accu, 
              precision = r, shifts = list(location = a, scale = b), z = atanh(rhoc), 
              var.z = var.z, confint = ci, level = level, center = xbar, cov = xcov, 
              bland = bland, equal.means = equal.means)
  if (ustat)
    out$ustat <- list(rhoc = u$rhoc, var.rhoc = u$var.rhoc, ustat = u$u, cov = u$v)
  if (equal.means) {
    out$Restricted <- list(ccc = rho0, var.ccc = var.rho0, accuracy = acc0, precision = r0, 
    shifts = list(location = 0, scale = b0), L1 = list(rho1 = 1 - sqrt(1 - rho0), 
    var.rho1 = 0.25 * var.rho0 / (1 - rho0)), center = rep(lambda, p), cov = Sigma)
  }
  class(out) <- "ccc"
  out
}

print.ccc <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  dput(x$call, control = NULL)
  cf <- c(x$ccc, x$var.ccc, x$accuracy, x$precision)
  names(cf) <- c("estimate","variance","accuracy","precision")
  cat("\nCoefficients:\n ")
  print(format(round(cf, digits = digits)), quote = F, ...)
  cat("\n")
  cat(paste("Asymptotic ", 100 * x$level, "% confidence interval:\n", sep = ""))
  print(format(round(x$confint, digits = digits)), quote = F, ...)
  if (x$equal.means) {
    c0 <- c(x$Restricted$ccc, x$Restricted$var.ccc, x$Restricted$accuracy, x$Restricted$precision)
    names(c0) <- c("estimate","variance","accuracy","precision")
    cat("\nCoefficients under equality of means:\n ")
    print(format(round(c0, digits = digits)), quote = F, ...)
  }
  invisible(x)
}

ustat.rhoc <- function(x)
{ ## estimation and asymptotic variance of rhoc based on U-statistics
  z <- x
  x <- z[,1]
  y <- z[,2]
  n <- nrow(z)

  # computing underlying kernels for the U-statistics
  z <- .Fortran("rhoc_ustat",
                x = as.double(x),
                y = as.double(y),
                n = as.integer(n),
                p1 = double(n),
                p2 = double(n),
                p3 = double(n))

  # estimating U-statistics and their covariance
  phi <- cbind(z$p1, z$p2, z$p3)
  ones <- rep(1, n)
  p <- 3
  storage.mode(phi) <- "double"
  z <- .C("cov_weighted",
          phi = phi,
          n = as.integer(n),
          p = as.integer(p),
          weights = as.double(ones),
          mean = double(p),
          cov  = double(p * p))[c("mean","cov")]
  u <- z$mean
  v <- (4 / n) * matrix(z$cov, nrow = p, ncol = p)

  # estimating rhoc using U-statistics
  u1 <- u[1]
  u2 <- u[2]
  u3 <- u[3]
  h = (n - 1) * (u3 - u1)
  g = u1 + n * u2 + (n - 1) * u3
  rhoc = h / g

  # asymptotic variance of rhoc based on U-statistics
  v11 <- v[1,1]
  v12 <- v[1,2]
  v13 <- v[1,3]
  v22 <- v[2,2]
  v23 <- v[2,3]
  v33 <- v[3,3]
  var.h = (v11 + v33 - 2 * v13) * (n - 1)^2
  var.g = v33 * (n - 1)^2 + v11 + v22 * n^2 + 2 * (n - 1) * v13 + 2 * n * (n - 1) * v23 + 2 * n * v12
  cov.hg = (n - 1) * ((n - 1) * v33 + n * v23  - (n - 2) * v13 - v11 - n  * v12)
  var.rhoc = rhoc^2 * (var.h / h^2 + var.g / g^2 - 2 * cov.hg / (h * g))

  # output object  
  o <- list(rhoc = rhoc, var.rhoc = var.rhoc, kernel = phi, u = u, v = v)
  o
}
