## ID: ols.R, last updated 2022-10-06, F.Osorio

ols <-
function(formula, data, subset, na.action, method = "qr", tol = 1e-7, maxiter = 100,
  model = FALSE, x = FALSE, y = FALSE, contrasts = NULL, ...)
{ ## ordinary least-squares fit
  ret.x <- x
  ret.y <- y
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$method <- mf$tol <- mf$maxiter <- mf$model <- mf$x <- mf$y <- NULL
  mf$contrasts <- mf$... <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(Terms, mf, contrasts)

  # call fitter
  z <- ols.fit(x, y, method, tol, maxiter)

  # output
  z$call <- Call
  z$na.action <- attr(mf, "na.action")
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(Terms, mf)
  z$terms <- Terms
  if (model)
    z$model <- mf
  if (ret.y)
    z$y <- y
  if (ret.x)
    z$x <- x
  class(z) <- "ols"
  z
}

ols.fit <-
function(x, y, method = "qr", tol = 1e-7, maxiter = 100)
{ ## dispatcher among various fitting functions
  if (!is.numeric(x))
    stop("model matrix must be numeric.")
  if (!is.numeric(y))
    stop("response must be numeric.")
  if (!length(x))
    method <- "null"
  contr <- attr(x, 'contrast')
  fit <- switch(method,
                cg    = ols.fit.cg(x, y, tol, maxiter),
                chol  = ols.fit.chol(x, y),
                qr    = ols.fit.qr(x, y),
                svd   = ols.fit.svd(x, y),
                sweep = ols.fit.sweep(x, y),
                stop(paste("unimplemented method:", method)))
  fit$contrast <- contr
  fit$method <- method
  fit
}

ols.fit.chol <-
function(x, y)
{
  if (is.matrix(y))
    stop("'ols.fit.chol' does not support multiple responses")
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
	ny <- length(y)
	dx <- dim(x)
	n <- dx[1]
	xn <- dimnames(x)[[2]]
	if (n != ny)
		stop("Number of observations in x and y not equal")
	p <- dx[2]
	xx <- crossprod(x)
	xy <- crossprod(x, y)

  z <- .C("chol_dcmp",
          R = xx,
          ldR = as.integer(p),
          p = as.integer(p),
          job = as.integer(1),
          info = as.integer(0))[c("R","info")]
  if (z$info) {
    if (z$info > 0)
      stop("Cholesky decomposition not of full rank")
  }

  R <- z$R # is not upper triangular
  eff <- backsolve(R, xy, transpose = TRUE)
  b <- backsolve(R, eff)
  fitted <- x %*% b
	fitted <- fitted[,]
  residuals <- y - fitted
	z <- list(coefficients = b[,], residuals = residuals, fitted.values = fitted,
    RSS = minkowski(residuals)^2, cov.unscaled = chol2inv(R), dims = dx)
  names(z$coefficients) <- xn
  dimnames(z$cov.unscaled) <- list(xn, xn)
  class(z) <- "ols"
	z
}

ols.fit.cg <-
function(x, y, tol = 1e-7, maxiter = 100)
{
  if (is.matrix(y))
    stop("'ols.fit.gc' does not support multiple responses")
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
	ny <- length(y)
	dx <- dim(x)
	n <- dx[1]
	xn <- dimnames(x)[[2]]
	if (n != ny)
		stop("Number of observations in x and y not equal")
	p <- dx[2]
  now <- proc.time()

  z <- .C("OLS_cg",
          x = x,
          ldx = as.integer(n),
          n = as.integer(n),
          p = as.integer(p),
          y = y,
          coef = double(p),
          tol = as.double(tol),
          maxiter = as.integer(maxiter),
          info = as.integer(0))[c("coef","info")]

  speed <- proc.time() - now
  fitted <- x %*% z$coef
	fitted <- fitted[,]
  residuals <- y - fitted
	z <- list(coefficients = z$coef, residuals = residuals, fitted.values = fitted,
    RSS = minkowski(residuals)^2, cov.unscaled = NULL, dims = dx, iter = z$info,
    speed = speed)
  names(z$coefficients) <- xn
  class(z) <- "ols"
	z
}

ols.fit.qr <-
function(x, y)
{
  if (is.matrix(y))
    stop("'ols.fit.qr' does not support multiple responses")
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
	ny <- length(y)
	dx <- dim(x)
	n <- dx[1]
	xn <- dimnames(x)[[2]]
	if (n != ny)
		stop("Number of observations in x and y not equal")
	p <- dx[2]

  z <- .C("OLS_qr",
          x = x,
          ldx = as.integer(n),
          n = as.integer(n),
          p = as.integer(p),
          y = y,
          qraux = double(p),
          coef = double(p),
          fitted = double(n),
          resid = double(n),
          RSS = as.double(0))[c("x","qraux","coef","fitted","resid","RSS")]
  R <- z$x[1:p,1:p] # is not upper triangular
	z <- list(coefficients = z$coef, residuals = z$resid, fitted.values = z$fitted,
    RSS = z$RSS, cov.unscaled = chol2inv(R), dims = dx)
  names(z$coefficients) <- xn
  names(z$residuals) <- 1:n
  names(z$fitted.values) <- 1:n
  dimnames(z$cov.unscaled) <- list(xn, xn)
  class(z) <- "ols"
	z
}

ols.fit.svd <-
function(x, y)
{
  if (is.matrix(y))
    stop("'ols.fit.svd' does not support multiple responses")
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
	ny <- length(y)
	dx <- dim(x)
  xn <- dimnames(x)[[2]]
  n <- dx[1]
  p <- dx[2]
  if (n < p)
    stop("'ols.fit.svd' only implemented for overdetermined systems")

  v <- double(p * p)
  dim(v) <- c(p, p)
  job <- 21 # using the storage of 'x'
  z <- .C("svd_dcmp",
          x = x,
          ldx = as.integer(n),
          n = as.integer(n),
          p = as.integer(p),
          u = x,
          ldu = as.integer(n),
          d = double(p),
          v = v,
          ldv = as.integer(p),
          job = as.integer(job),
          info = as.integer(0))[c("u","d","v","info")]
	if (z$info)
		stop(paste("Numerical error (code", z$info, ") in algorithm"))

	u <- crossprod(z$u, y)
  u <- u[,]
	d <- z$d
	if (all(dp <- d > 0))
		inv <- 1.0 / d
	else {
    inv <- d
		inv[dp] <- 1.0 / d[dp]
		inv[!dp] <- 0
	}
  # coefficients
  a <- inv * u[1:p]
	b <- z$v %*% a
  fitted <- x %*% b
	fitted <- fitted[,]
  residuals <- y - fitted
  z <- list(coefficients = b[,], residuals = residuals, fitted.values = fitted,
    RSS = minkowski(residuals)^2, cov.unscaled = z$v %*% (inv^2 * t(z$v)),
    dims = dx, rank = sum(dp))
  names(z$coefficients) <- xn
  dimnames(z$cov.unscaled) <- list(xn, xn)
  class(z) <- "ols"
	z
}

ols.fit.sweep <-
function(x, y)
{
  if (is.matrix(y))
    stop("'ols.fit.sweep' does not support multiple responses")
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
	ny <- length(y)
	dx <- dim(x)
	n <- dx[1]
  xn <- dimnames(x)[[2]]
	if (n != ny)
		stop("Number of observations in x and y not equal")
	p <- dx[2]
	z <- crossprod(cbind(x, y))
  pm1 <- p + 1

  z <- .C("sweep_operator",
          z = z,
          ldz = as.integer(pm1),
          pm1 = as.integer(pm1),
          k = as.integer(1:p),
          p = as.integer(p),
          reverse = as.integer(0))$z

  b <- z[1:p,pm1]
  fitted <- x %*% b
	fitted <- fitted[,]
  z <- list(coefficients = b, residuals = y - fitted, fitted.values = fitted,
    RSS = z[pm1,pm1], cov.unscaled = z[1:p,1:p], dims = dx)
  names(z$coefficients) <- xn
  dimnames(z$cov.unscaled) <- list(xn, xn)
	class(z) <- "ols"
	z
}

## extractors
deviance.ols <- function(object, ...) object$RSS
residuals.ols <- function(object, ...) object$residuals

print.ols <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)
  nobs <- x$dims[1]
  rdf <- nobs - x$dims[2]
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
  if (rdf > 0) {
    RSS <- x$RSS
    cat("Residual standard error:", format((RSS / rdf)^0.5), "\n")
  }
  invisible(x)
}

summary.ols <- function (object, correlation = FALSE, ...)
{
  z <- object
  n <- z$dims[1]
  p <- z$dims[2]
  rdf <- n - p
  if (p == 0) {
    r <- z$residuals
    RSS <- z$RSS
    resvar <- RSS / rdf
    ans <- z[c("call", "terms")]
    class(ans) <- "summary.ols"
    ans$residuals <- r
    ans$df <- c(0L, n)
    ans$coefficients <- matrix(NA_real_, 0L, 4L, dimnames =
      list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
    ans$sigma <- sqrt(resvar)
    ans$cov.unscaled <- matrix(NA_real_, 0L, 0L)
    if (correlation)
      ans$correlation <- ans$cov.unscaled
    return(ans)
  }
  if (is.null(z$terms))
    stop("invalid 'ols' object:  no 'terms' component")

  r <- z$residuals
  f <- z$fitted.values
  RSS <- z$RSS
  resvar <- RSS / rdf
  if (is.finite(resvar) && resvar < (mean(f)^2 + var(c(f))) * 1e-30) # a few times .Machine$double.eps^2
    warning("essentially perfect fit: summary may be unreliable")
  ans <- z[c("call", "terms")]
  if (z$method != "cg") {
    se <- sqrt(diag(z$cov.unscaled) * resvar)
    est <- z$coefficients
    tval <- est / se
    ans$coefficients <- cbind(Estimate = est, "Std. Error" = se, "t value" = tval,
                              "Pr(>|t|)" = 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    ans$cov.unscaled <- z$cov.unscaled
    if (correlation)
      ans$correlation <- (ans$cov.unscaled * resvar) / outer(se, se)
  } else {
    ans$coefficients <- z$coefficients
    ans$cov.unscaled <- NULL
    if (correlation)
      ans$correlation <- NULL
  }
  ans$residuals <- r
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf)
  ans$logLik <- logLik(z)
  class(ans) <- "summary.ols"
  ans
}

print.summary.ols <-
function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep = "")
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  cat("Residuals:\n")
  if (rdf > 5L) {
    rq <- quantile(resid)
    names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(rq, digits = digits, ...)
  } else if (rdf > 0L) {
    print(resid, digits = digits, ...)
  } else { # rdf == 0 : perfect fit!
    cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
    cat("\n")
  }
  cat("\nCoefficients:\n")
  print(format(round(x$coef, digits = digits)), quote = FALSE, ...)
  cat("\nResidual standard error:", format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom")
  cat("\n")
	cat("Log-likelikood:", format(signif(x$logLik, digits)))
  cat("\n")
  invisible(x)
}

logLik.ols <- function(object, ...)
{ # log-likelihood for ols objects
  res <- object$residuals # not resid(object) because of NA methods
  n <- object$dims[1]
  p <- object$dims[2]
  RSS <- object$RSS
  val <- -.5 * n * (log(2 * pi) + 1 - log(n) + log(RSS))
  attr(val, "nobs") <- n
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
}
