\name{symm.prod}
\alias{symm.prod}
\title{Matrix multiplication envolving the symmetrizer matrix}
\description{
  Given the order of a symmetrizer matrix \eqn{\bold{N}} of order \eqn{n} and a
  conformable matrix \eqn{\bold{X}}, performs one of the matrix-matrix operations:
  \itemize{
    \item \eqn{\bold{Y} = \bold{NX}}, if \code{side = "left"}, or
    \item \eqn{\bold{Y} = \bold{XN}}, if \code{side = "right"},
  }
  The main aim of \code{symm.prod} is to do this matrix multiplication \bold{without forming}
  the symmetrizer matrix.
}
\usage{
symm.prod(n = 1, x = NULL, side = "left")
}
\arguments{
  \item{n}{ order of the symmetrizer matrix.}
  \item{x}{ numeric matrix (or vector).}
  \item{side}{ a string selecting if symmetrizer matrix is pre-multiplying \eqn{\bold{X}},
  that is \code{side = "left"} or post-multiplying \eqn{\bold{X}}, by using \code{side = "right"}.}
}
\details{
  Underlying \code{C} code only uses information provided by \code{\link{symm.info}} to
  performs the matrix multiplication. The symmetrizer matrix is \bold{never} created.
}
\seealso{
  \code{\link{symmetrizer}}
}
\examples{
N4 <- symmetrizer(n = 4, matrix = TRUE)
x <- matrix(1:32, ncol = 2)
y <- N4 \%*\% x

z <- symm.prod(n = 4, x) # N4 is not stored
all(z == y) # matrices y and z are equal!
}
\keyword{array}
\keyword{algebra}
