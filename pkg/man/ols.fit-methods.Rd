\name{ols.fit-methods}
\alias{ols.fit.cg}
\alias{ols.fit.chol}
\alias{ols.fit.qr}
\alias{ols.fit.svd}
\alias{ols.fit.sweep}
\title{Fit a linear model}
\usage{
ols.fit.cg(x, y, tol = 1e-7, maxiter = 100)
ols.fit.chol(x, y)
ols.fit.qr(x, y)
ols.fit.svd(x, y)
ols.fit.sweep(x, y)
}
\description{
  Fits a linear model, returning the bare minimum computations.
}
\arguments{
  \item{x, y}{ numeric vectors or matrices for the predictors and the response in
  a linear model. Typically, but not necessarily, \code{x} will be constructed by
  one of the fitting functions.}
  \item{tol}{ tolerance for the conjugate gradients (\code{gc}) method. Default is
  \code{tol = 1e-7}.}
  \item{maxiter}{ The maximum number of iterations for the conjugate gradients (\code{gc})
  method. Defaults to 100.}
}
\value{
  The bare bones of an \code{ols} object: the coefficients, residuals, fitted values,
  and some information used by \code{summary.ols}.
}
\seealso{
  \code{\link{ols}}, \code{\link{ols.fit}}, \code{\link{lm}}
}
\examples{
set.seed(151)
n <- 100
p <- 2
x <- matrix(rnorm(n * p), n, p) # no intercept!
y <- rnorm(n)
z <- ols.fit.chol(x, y)
z
}
\keyword{regression}
\keyword{array}
