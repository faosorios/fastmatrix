\name{whitening}
\alias{whitening}
\title{Whitening transformation}
\description{
  Applies the whitening transformation to a data matrix based on the Cholesky decomposition
  of the empirical covariance matrix.
}
\usage{
whitening(x, Scatter = NULL)
}
\arguments{
  \item{x}{ vector or matrix of data with, say, \eqn{p} columns.}
  \item{Scatter}{ covariance (or scatter) matrix (\eqn{p \times p}{p x p}) of the
    distribution, must be positive definite. If \code{NULL}, the covariance matrix
    is estimated from the data.}
}
\value{
  Returns the whitened data matrix \eqn{\bold{Z} = \bold{X W}^T}, where
  \deqn{\bold{W}^T\bold{W} = \bold{S}^{-1},}
  with \eqn{\bold{S}} the empirical covariance matrix.
}
\references{
  Kessy, A., Lewin, A., Strimmer, K. (2018).
  Optimal whitening and decorrelation.
  \emph{The American Statistician} \bold{72}, 309-314.
}
\examples{
x <- iris[,1:4]
species <- iris[,5]
pairs(x, col = species) # plot of Iris

# whitened data
z <- whitening(x)
pairs(z, col = species) # plot of
}
\keyword{multivariate}
