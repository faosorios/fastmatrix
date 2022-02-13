\name{cholupdate}
\alias{cholupdate}
\title{Rank 1 update to Cholesky factorization}
\usage{
cholupdate(r, x)
}
\arguments{
  \item{r}{ a upper triangular matrix, the Cholesky factor of matrix a.}
  \item{x}{ vector defining the rank one update.}
}
\description{
  function \code{chol.update}, where \code{R = chol(A)} is the original Cholesky
  factorization of \eqn{\bold{A}}, returns the upper triangular Cholesky factor of
  \eqn{\bold{A} + \bold{xx}^T}, with \code{x} a column vector of appropriate dimension.
}
\references{
  Golub, G.H., Van Loan, C.F. (2013).
  \emph{Matrix Computations}, 4th Edition.
  John Hopkins University Press.
}
\seealso{
  \code{\link{chol}}
}
\examples{
a <- matrix(c(1,1,1,1,2,3,1,3,6), ncol = 3)
r <- chol(a)
x <- c(0,0,1)
b <- a + outer(x,x)
r1 <- cholupdate(r, x)
r1
all(r1 == chol(b)) # TRUE
}
\keyword{algebra}
\keyword{array}
