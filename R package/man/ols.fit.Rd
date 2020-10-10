\name{ols.fit}
\alias{ols.fit}
\title{Fitter Functions for Linear Models}
\usage{
ols.fit(x, y, method = "qr", \dots)
}
\description{
  This function is a \emph{switcher} among various numerical fitting functions
  (\code{\link{ols.fit.chol}}, \code{\link{ols.fit.qr}}, \code{\link{ols.fit.svd}}
  and \code{\link{ols.fit.sweep}}). The argument \code{method} does the switching:
  \code{"qr"} for \code{\link{ols.fit.qr}}, etc. This should usually \emph{not} be
  used directly unless by experienced users.
}
\arguments{
  \item{x}{ design matrix of dimension \code{n * p}.}
  \item{y}{ vector of observations of length \code{n}.}
  \item{method}{ currently, methods \code{"chol"}, \code{"qr"} (default), \code{"svd"}
  and \code{"sweep"} are supported.}
  \item{\dots}{ currently disregarded.}
}
\value{
  a \code{\link{list}} with components:
  \item{coefficients}{ a named vector of coefficients}
  \item{residuals}{ the residuals, that is response minus fitted values.}
  \item{fitted.values}{ the fitted mean values.}
  \item{RSS}{ the residual sum of squares.}
  \item{cov.unscaled}{ a \eqn{p \times p}{p x p} matrix of (unscaled)
    covariances of the \eqn{\hat\beta_j}{coef[j]}, \eqn{j=1, \dots, p}.}
}
\seealso{
  \code{\link{ols.fit.chol}}, \code{\link{ols.fit.qr}}, \code{\link{ols.fit.svd}},
  and \code{\link{ols.fit.sweep}}.
}
\examples{
set.seed(151)
n <- 100
p <- 2
x <- matrix(rnorm(n * p), n, p) # no intercept!
y <- rnorm(n)
fm <- ols.fit(x = x, y = y, method = "chol")
fm
}
\keyword{regression}
\keyword{array}
