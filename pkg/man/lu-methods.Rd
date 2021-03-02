\name{lu-methods}
\alias{constructX}
\alias{extractL}
\alias{extractU}
\title{Reconstruct the L, U, or X Matrices from an LU object}
\usage{
constructX(x)
extractL(x)
extractU(x)
}
\arguments{
  \item{x}{object representing an LU factorization. This will typically have
  come from a previous call to \code{\link{lu}}.}
}
\description{
  Returns the original matrix from which the object was constructed or
  the components of the factorization.
}
\value{
  \code{constructX} returns \eqn{\bold{X}}, the original matrix from which the \code{lu}
  object was constructed (because of the pivoting the \eqn{\bold{X}} matrix is not exactly
  the product between \eqn{\bold{L}} and \eqn{\bold{U}}).

  \code{extractL} returns \eqn{\bold{L}}. This may be pivoted.

  \code{extractU} returns \eqn{\bold{U}}.
}
\seealso{
  \code{\link{lu}}.
}
\examples{
a <- matrix(c(10,-3,5,-7,2,-1,0,6,5), ncol = 3)
z <- lu(a)
L <- extractL(z)
L
U <- extractU(z)
U
X <- constructX(z)
all(a == X)
}
\keyword{algebra}
\keyword{array}
