\name{ols}
\alias{ols}
\title{Fit linear regression model}
\description{
  Returns an object of class \code{"ols"} that represents a linear model fit.
}
\usage{
ols(formula, data, subset, na.action, method = "qr", tol = 1e-7, maxiter = 100,
  x = FALSE, y = FALSE, contrasts = NULL, ...)
}
\arguments{
  \item{formula}{ an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the model to
    be fitted.}
  \item{data}{ an optional data frame, list or environment (or object coercible
    by \code{\link{as.data.frame}} to a data frame) containing the variables in
    the model. If not found in \code{data}, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{ols}
    is called.}
  \item{subset}{ an optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{na.action}{ a function which indicates what should happen when the data
    contain \code{NA}s. The default is set by the \code{na.action} setting of
    \code{\link{options}}, and is \code{\link{na.fail}} if that is unset.}
  \item{method}{the least squares fitting method to be used; the options are \code{"cg"}
    (conjugate gradients), \code{"chol"}, \code{"qr"} (the default), \code{"svd"} and
    \code{"sweep"}.}
  \item{tol}{ tolerance for the conjugate gradients (\code{gc}) method. Default is
    \code{tol = 1e-7}.}
  \item{maxiter}{ The maximum number of iterations for the conjugate gradients (\code{gc})
    method. Defaults to 100.}
  \item{x, y}{ logicals.  If \code{TRUE} the corresponding components of
    the fit (the model matrix, the response) are returned.}
  \item{contrasts}{ an optional list. See the \code{contrasts.arg} of
    \code{\link{model.matrix.default}}.}
  \item{\dots}{ additional arguments (currently disregarded).}
}
\value{
  \code{ols} returns an object of \code{\link{class}} \code{"ols"}.

  The function \code{summary} is used to obtain and print a summary of the
  results.  The generic accessor functions \code{coefficients}, \code{fitted.values}
  and \code{residuals} extract various useful features of the value returned by \code{ols}.

  An object of class \code{"ols"} is a list containing at least the
  following components:

  \item{coefficients}{a named vector of coefficients}
  \item{residuals}{the residuals, that is response minus fitted values.}
  \item{fitted.values}{the fitted mean values.}
  \item{RSS}{the residual sum of squares.}
  \item{cov.unscaled}{a \eqn{p \times p}{p x p} matrix of (unscaled)
    covariances of the \eqn{\hat\beta_j}{coef[j]}, \eqn{j=1, \dots, p}.}
  \item{call}{the matched call.}
  \item{terms}{the \code{\link{terms}} object used.}
  \item{contrasts}{(only where relevant) the contrasts used.}
  \item{y}{if requested, the response used.}
  \item{x}{if requested, the model matrix used.}
  \item{model}{if requested (the default), the model frame used.}
}
\seealso{
  \code{\link{ols.fit}}, \code{\link{lm}}, \code{\link{lsfit}}
}
\examples{
# tiny example of regression
y <- c(1, 3, 3, 2, 2, 1)
x <- matrix(c(1, 1,
              2, 1,
              3, 1,
              1,-1,
              2,-1,
              3,-1), ncol = 2, byrow = TRUE)
f0 <- ols(y ~ x) # intercept is included by default
f0 # printing results (QR method was used)

f1 <- ols(y ~ x, method = "svd") # using SVD method instead
f1
}
\keyword{regression}
