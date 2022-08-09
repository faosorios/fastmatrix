\name{vech}
\alias{vech}
\title{Vectorization the lower triangular part of a square matrix}
\description{
  This function returns a vector obtained by stacking the lower triangular part of a
  square matrix.
}
\usage{
vech(x)
}
\arguments{
  \item{x}{ a square matrix.}
}
\value{
  Let \eqn{\bold{X}} be a \eqn{n} by \eqn{n} matrix, then \code{vech}(\eqn{\bold{X}})
  is a \eqn{n(n+1)/2}-dimensional vector.
}
\examples{
x <- matrix(rep(1:10, each = 10), ncol = 10)
x
y <- vech(x)
y
}
\keyword{array}
