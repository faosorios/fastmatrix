\name{ridge}
\alias{ridge}
\title{Ridge regression}
\description{
  Fit a linear model by ridge regression, returning an object of class \code{"ridge"}.
}
\usage{
ridge(formula, data, subset, lambda = 1.0, method = "GCV", ngrid = 200, tol = 1e-07,
  maxiter = 50, na.action, model = FALSE, x = FALSE, y = FALSE, contrasts = NULL, ...)
}
\arguments{
  \item{formula}{ an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the model to
    be fitted.}
  \item{data}{ an optional data frame, list or environment (or object coercible
    by \code{\link{as.data.frame}} to a data frame) containing the variables in
    the model. If not found in \code{data}, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{ridge}
    is called.}
  \item{subset}{ an optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{na.action}{ a function which indicates what should happen when the data
    contain \code{NA}s. The default is set by the \code{na.action} setting of
    \code{\link{options}}, and is \code{\link{na.fail}} if that is unset.}
  \item{lambda}{ a scalar or vector of ridge constants. A value of 0 corresponds
  to ordinary least squares.}
  \item{method}{ the method for choosing the ridge parameter lambda. If \code{method = "none"},
  then lambda is 'fixed'. If \code{method = "GCV"} (the default) then the ridge
  parameter is chosen automatically using the generalized cross validation (GCV) criterion.
  For \code{method = "grid"}, optimal value of lambda is selected computing the GCV
  criterion over a grid. If \code{method = "MSE"} the optimal ridge parameter is selected
  minimizing the mean squared estimation error criterion, this is the \code{ORPS1}
  subroutine by Lee (1987).}
  \item{ngrid}{number of elements in the grid used to compute the GCV criterion.
  Only required if \code{method = "grid"} and \code{lambda} is a scalar.}
  \item{tol}{ tolerance for the optimization of the GCV criterion. Default is \code{1e-7}.}
  \item{maxiter}{ maximum number of iterations. The default is 50.}
  \item{model, x, y}{ logicals.  If \code{TRUE} the corresponding components of
    the fit (the model frame, the model matrix, the response) are returned.}
  \item{contrasts}{ an optional list. See the \code{contrasts.arg} of
    \code{\link{model.matrix.default}}.}
  \item{\dots}{ additional arguments to be passed to the low level regression
    fitting functions (not implemented).}
}
\details{
  \code{ridge} function fits in linear ridge regression \bold{without} scaling or centering
  the regressors and the response. In addition, If an intercept is present in the model, its
  coefficient is penalized.)
}
\value{
  A list with the following components:
  \item{dims}{ dimensions of model matrix.}
  \item{coefficients}{ a named vector of coefficients.}
  \item{scale}{ a named vector of coefficients.}
  \item{fitted.values}{ the fitted mean values.}
  \item{residuals}{ the residuals, that is response minus fitted values.}
  \item{RSS}{ the residual sum of squares.}
  \item{edf}{ the effective number of parameters.}
  \item{GCV}{ vector (if \code{method = "grid"}) of GCV values.}
  \item{HKB}{ HKB estimate of the ridge constant.}
  \item{LW}{ LW estimate of the ridge constant.}
  \item{lambda}{ vector (if \code{method = "grid"}) of lambda values; otherwise, for
  methods \code{method = "none"}, \code{"GCV"} or \code{"MSE"}, the value of ridge
  parameter used by the algorithm.}
  \item{optimal}{ value of lambda with the minimum GCV (only relevant if \code{method = "grid"}).}
  \item{iterations}{ number of iterations performed by the algorithm (only relevant if \code{method = "MSE"}).}
  \item{call}{ the matched call.}
  \item{terms}{ the \code{\link{terms}} object used.}
  \item{contrasts}{ (only where relevant) the contrasts used.}
  \item{y}{ if requested, the response used.}
  \item{x}{ if requested, the model matrix used.}
  \item{model}{ if requested, the model frame used.}
}
\references{
  Golub, G.H., Heath, M., Wahba, G. (1979).
  Generalized cross-validation as a method for choosing a good ridge parameter.
  \emph{Technometrics} \bold{21}, 215-223.

  Hoerl, A.E., Kennard, R.W., Baldwin, K.F. (1975).
  Ridge regression: Some simulations.
  \emph{Communication in Statistics} \bold{4}, 105-123.

  Hoerl, A.E., Kennard, R.W. (1970).
  Ridge regression: Biased estimation of nonorthogonal problems.
  \emph{Technometrics} \bold{12}, 55-67.

  Lawless, J.F., Wang, P. (1976).
  A simulation study of ridge and other regression estimators.
  \emph{Communications in Statistics} \bold{5}, 307-323.

  Lee, T.S (1987).
  Algorithm AS 223: Optimum ridge parameter selection.
  \emph{Applied Statistics} \bold{36}, 112-118.
}
\seealso{
  \code{\link{lm}}, \code{\link{ols}}
}
\examples{
z <- ridge(GNP.deflator ~ ., data = longley, lambda = 4, method = "grid")
z # ridge regression on a grid over seq(0, 4, length = 200)

z <- ridge(GNP.deflator ~ ., data = longley)
z # ridge parameter selected using GCV (default)
}
\keyword{models}
