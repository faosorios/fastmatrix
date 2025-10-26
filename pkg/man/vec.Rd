\name{vec}
\alias{vec}
\title{Vectorization of a matrix}
\description{
  This function returns a vector obtained by stacking the columns of \eqn{\bold{X}}.
}
\usage{
vec(x)
}
\arguments{
  \item{x}{ a numeric matrix.}
}
\value{
  Let \eqn{\bold{X}} be a \eqn{n} by \eqn{m} matrix, then \code{vec}(\eqn{\bold{X}})
  is a \eqn{nm}-dimensional vector.
}
\examples{
x <- matrix(rep(1:10, each = 10), ncol = 10)
x
y <- vec(x)
y
}
\keyword{array}
