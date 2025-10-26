\name{asSymmetric}
\alias{asSymmetric}
\title{Force a matrix to be symmetric}
\description{
  Force a square matrix \eqn{\bold{X}} to be symmetric
}
\usage{
asSymmetric(x, lower = TRUE)
}
\arguments{
  \item{x}{ a square matrix to be forced to be symmetric.}
  \item{lower}{ logical, should the upper (lower) triangle be replaced with the
  lower (upper) triangle?}
}
\value{
  a square symmetric matrix.
}
\examples{
a <- matrix(1:16, ncol = 4)
isSymmetric(a) # FALSE
a <- asSymmetric(a) # copy lower triangle into upper triangle
}
\keyword{array}
