\name{Mahalanobis}
\alias{Mahalanobis}
\title{Mahalanobis distance}
\usage{
Mahalanobis(x, center, cov, inverted = FALSE)
}
\arguments{
  \item{x}{vector or matrix of data. As usual, rows are observations and columns are
    variables.}
  \item{center}{mean vector of the distribution.}
  \item{cov}{covariance matrix (\eqn{p \times p}{p x p}) of the distribution, must
    be positive definite.}
  \item{inverted}{logical. If \code{TRUE}, \code{cov} is supposed to contain the
    \emph{inverse} of the covariance matrix.}
}
\description{
  Returns the squared Mahalanobis distance of all rows in \eqn{\bold{x}} and the
  vector \eqn{\bold{\mu}} = \code{center} with respect to \eqn{\bold{\Sigma}} = \code{cov}.
  This is (for vector \eqn{\bold{x}}) defined as
  \deqn{D^2 = (\bold{x} - \bold{\mu})^T \bold{\Sigma}^{-1} (\bold{x} - \bold{\mu})}
}
\details{
  Unlike function \code{mahalanobis}, the covariance matrix is factorized using the
  Cholesky decomposition, which allows to assess if \code{cov} is positive definite.
  Unsuccessful results from the underlying LAPACK code will result in an error message.
}
\seealso{\code{\link{cov}}, \code{\link{mahalanobis}}}
\examples{
x <- cbind(1:6, 1:3)
xbar <- colMeans(x)
S <- matrix(c(1,4,4,1), ncol = 2) # is negative definite
D2 <- mahalanobis(x, center = xbar, S)
all(D2 >= 0) # several distances are negative

## next command produces the following error:
## Covariance matrix is possibly not positive-definite
\dontrun{D2 <- Mahalanobis(x, center = xbar, S)}
}
\keyword{multivariate}
