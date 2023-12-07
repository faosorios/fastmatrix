\name{mchol}
\alias{mchol}
\title{The modified Cholesky factorization}
\usage{
mchol(x)
}
\arguments{
  \item{x}{ a symmetric but not necessarily positive definite matrix to be factored.}
}
\description{
  Compute the Cholesky factorization of a real symmetric but not necessarily positive 
  definite matrix.
}
\value{
  The lower triangular factor of modified Cholesky decomposition, i.e., the matrix 
  \eqn{\bold{L}} such that \eqn{\bold{X} + \bold{E} = \bold{LL}^T}, where \eqn{\bold{E}} 
  is a nonnegative diagonal matrix that is zero if \eqn{\bold{X}} es sufficiently positive 
  definite.
}
\references{
  Gill, P.E., Murray, W., Wright, M.H. (1981).
  \emph{Practical Optimization}.
  Academic Press, London.

  Nocedal, J., Wright, S.J. (1999).
  \emph{Numerical Optimization}.
  Springer, New York.
}
\seealso{
  \code{\link{chol}}, \code{\link{ldl}}
}
\examples{
# a non-positive-definite matrix
a <- matrix(c(4,2,1,2,6,3,1,3,-.004), ncol = 3)
try(chol(a)) # fails
z <- mchol(a)
z # triangular factor

# modified 'a' matrix
tcrossprod(z)
}
\keyword{algebra}
\keyword{array}
