\name{vec}
\alias{vec}
\title{vectorization of a matrix}
\description{
  This function returns a vector obtained by stacking the columns of x
}
\usage{
  vec(x)
}
\arguments{
  \item{x}{ a numeric matrix.}
}
\value{
  Let \code{x} be a \eqn{n} by \eqn{m} matrix, then \code{vec(x)} is a \eqn{nm}-dimensional
  vector.
}
\examples{
x <- matrix(rep(1:10, each = 10), ncol = 10)
x
y <- vec(x)
y
}
\keyword{array}
