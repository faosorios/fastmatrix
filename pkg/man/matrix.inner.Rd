\name{matrix.inner}
\alias{matrix.inner}
\title{Compute the inner product between two rectangular matrices}
\description{
  Computes the inner product between two rectangular matrices calling BLAS.
}
\usage{
matrix.inner(x, y = x)
}
\arguments{
  \item{x}{ a numeric matrix.}
  \item{y}{ a numeric matrix.}
}
\value{
  a real value, indicating the inner product between two matrices.
}
\examples{
x <- matrix(c(1, 1, 1,
              1, 2, 1,
              1, 3, 1,
              1, 1,-1,
              1, 2,-1,
              1, 3,-1), ncol = 3, byrow = TRUE)
y <- matrix(1, nrow = 6, ncol = 3)
matrix.inner(x, y)

# must be equal
matrix.norm(x, type = "Frobenius")^2
matrix.inner(x)
}
\keyword{array}
\keyword{math}
