## ID: ridge.R, last updated 2022-02-07, F.Osorio

ridge <-
function(formula, data, subset, lambda = 1.0, method = "GCV", ngrid = 200, tol = 1e-07,
  maxiter = 50, na.action, model = FALSE, x = FALSE, y = FALSE, contrasts = NULL, ...)
{ ## ordinary ridge regression
  ret.x <- x
  ret.y <- y
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$lambda <- mf$method <- mf$ngrid <- mf$tol <- mf$maxiter <- NULL
  mf$model <- mf$x <- mf$y <- mf$contrasts <- mf$... <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(Terms, mf, contrasts)
  xn <- dimnames(x)[[2]]
  dx <- dim(x)
  n <- dx[1]
  p <- dx[2]
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
  method <- pmatch(method, c("none", "grid", "GCV", "MSE"))

  grid <- length(lambda)
  default <- lambda[1]
  switch(method,
         "none" = {
           task <- 0
           lambda <- default
           ngrid <- 1
           gcv <- 0.0
         },
         "grid" = {
           task <- 1
           if (grid > 1) {
             ngrid <- grid
             gcv <- double(ngrid)
           } else {
             lambda <- seq(0, to = lambda, length = ngrid)
             gcv <- double(ngrid)
           }
         },
         "GCV" = {
           task <- 2
           lambda <- 1.0
           ngrid <- 1
           gcv <- 0.0
         },
         "MSE" = {
           task <- 3
           lambda <- default
           ngrid <- 1
           gcv <- 0.0
         })

  # call fitter
  z <- .C("OLS_ridge",
          x = x,
          ldx = as.integer(n),
          n = as.integer(n),
          p = as.integer(p),
          y = y,
          coef = double(p),
          scale = as.double(0),
          fitted = double(n),
          resid = double(n),
          RSS = as.double(0),
          edf = as.double(0),
          pen = as.double(0),
          gcv = as.double(gcv),
          hkb = as.double(0),
          lw = as.double(0),
          lambda = as.double(lambda),
          opt = as.double(0),
          ngrid = as.integer(ngrid),
          task = as.integer(task),
          tolerance = as.double(tol),
          maxiter = as.integer(maxiter))

  # output
  z <- list(dims = c(n,p), coefficients = z$coef, scale = z$scale, fitted.values = z$fitted,
            residuals = z$resid, RSS = z$RSS, edf = z$edf, pen = z$pen, GCV = z$gcv, HKB = z$hkb,
            LW = z$lw, lambda = z$lambda, optimal = z$opt, iterations = z$maxiter)
  names(z$coefficients) <- xn
  z$call <- Call
  z$method <- method
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
  class(z) <- "ridge"
  z
}

print.ridge <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)
  cat("\n")
  switch(x$method,
         "none" = {
           cat("Ridge parameter:", format(round(x$lambda, 4)), "\n")
         },
         "grid" = {
           cat("Optimal ridge parameter:", format(round(x$optimal, 4)), "\n")
         },
         "GCV" = {
           cat("Estimated ridge parameter:", format(round(x$lambda, 4)), "\n")
         },
         "MSE" = {
           cat("Optimal ridge parameter:", format(round(x$lambda, 4)))
           cat(" (converged in", x$iterations, "iterations)\n")
         },
         "MSP" = {
           cat("Optimal ridge parameter:", format(round(x$lambda, 4)))
           cat(" (converged in", x$iterations, "iterations)\n")
         })
  cat("\nNumber of observations:", x$dims[1], "\n")
  cat("Effective number of parameters:", format(round(x$edf, 4)), "\n")
  cat("Scale parameter estimate:", format(round(x$scale, 4)), "\n")
  invisible(x)
}
