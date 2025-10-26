\name{kurtosis}
\alias{kurtosis}
\alias{skewness}
\title{Mardia's multivariate skewness and kurtosis coefficients}
\usage{
kurtosis(x)

skewness(x)
}
\arguments{
 \item{x}{ matrix of data with, say, \eqn{p} columns.}
}
\description{
  Functions to compute measures of multivariate skewness \eqn{(b_{1p})} and kurtosis
  \eqn{(b_{2p})} proposed by Mardia (1970),
  \deqn{b_{1p} = \frac{1}{n^2}\sum\limits_{i=1}^n\sum\limits_{j=1}^n ((\bold{x}_i -
  \overline{\bold{x}})^T\bold{S}^{-1}(\bold{x}_j - \overline{\bold{x}}))^3,}
  and
  \deqn{b_{2p} = \frac{1}{n}\sum\limits_{i=1}^n ((\bold{x}_i - \overline{\bold{x}})^T
  \bold{S}^{-1}(\bold{x}_i - \overline{\bold{x}}))^2.}
}
\references{
  Mardia, K.V. (1970).
  Measures of multivariate skewness and kurtosis with applications.
  \emph{Biometrika} \bold{57}, 519-530.

  Mardia, K.V., Zemroch, P.J. (1975).
  Algorithm AS 84: Measures of multivariate skewness and kurtosis.
  \emph{Applied Statistics} \bold{24}, 262-265.
}
\examples{
setosa <- iris[1:50,1:4]
kurtosis(setosa)
skewness(setosa)
}
\keyword{multivariate}
