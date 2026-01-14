\name{house}
\alias{house}
\title{Householder vector/reflection}
\description{
  This function returns the Householder matrix (also called Householder reflection) 
  or the Householder vector, which are constructed based on a nonzero vector \eqn{\bold{x}}. 
}
\usage{
house(x, matrix = FALSE)
}
\arguments{
  \item{x}{ an \eqn{n}-dimensional numeric vector.}
  \item{matrix}{ a logical indicating whether the commutation matrix will be returned.}
}
\details{
  A Householder transformation is a rank-1 modification of the identity matrix which can 
  be used to zero out selected elements of a vector. A Householder matrix \eqn{\bold{P}}
  adopts the form
  \deqn{\bold{P} - \tau\bold{uu}^T,}
  where \eqn{\bold{u}} is called Householder vector, and \eqn{\tau = 2/(\bold{u}^T\bold{u}).}
}
\value{
  If \code{matrix = TRUE} an \eqn{n} by \eqn{n} matrix is returned, otherwise \code{house} 
  function return a list containing the following components:

  \item{u}{ the Householder vector.}
  \item{tau}{ the scalar \eqn{\tau} required to construct the Householder reflector.}
}
\references{
  Golub, G.H., van Loan, C.F. (1996).
  \emph{Matrix Computations}, 3rd Edition. 
  The Johns Hopkins University Press, Baltimore.
}
\examples{
x <- c(3,1,5,1)
z <- house(x)
z

mat <- house(x, matrix = TRUE)
v <- mat \%*\% x
v[1,1] == minkowski(x) # TRUE
}
\keyword{array}
\keyword{algebra}
