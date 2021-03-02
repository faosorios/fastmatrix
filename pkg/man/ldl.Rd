\name{ldl}
\alias{ldl}
\title{The LDL decomposition}
\usage{
ldl(x)
}
\arguments{
  \item{x}{ a symmetric numeric matrix whose LDL decomposition is to be computed.}
}
\description{
  Compute the LDL decomposition of a real symmetric matrix.
}
\value{
  The factorization has the form \eqn{\bold{X} = \bold{LDL}^T}, where \eqn{\bold{D}}
  is a diagonal matrix and \eqn{\bold{L}} is unitary lower triangular.

  The LDL decomposition of \code{x} is returned as a list with components:
  \item{lower}{ the unitary lower triangular factor \eqn{\bold{L}}.}
  \item{d}{ a vector containing the diagonal elements of \eqn{\bold{D}}.}
}
\references{
  Golub, G.H., Van Loan, C.F. (1996).
  \emph{Matrix Computations}, 3rd Edition.
  John Hopkins University Press.
}
\seealso{
  \code{\link{chol}}
}
\examples{
a <- matrix(c(2,-1,0,-1,2,-1,0,-1,1), ncol = 3)
z <- ldl(a)
z # information of LDL factorization

# computing det(a)
prod(z$d) # product of diagonal elements of D

# a non-positive-definite matrix
m <- matrix(c(5,-5,-5,3), ncol = 2)
try(chol(m)) # fails
ldl(m)
}
\keyword{algebra}
\keyword{array}
