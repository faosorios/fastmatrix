\name{comm.prod}
\alias{comm.prod}
\title{Matrix multiplication envolving the commutation matrix}
\description{
  Given the row and column dimension of a commutation and matrix \code{x}, performs
  one of the matrix-matrix operations:
  \itemize{
    \item \eqn{{\bf{Y}} = {\bf{KX}}}, if \code{side = "left"} and \code{transposed = FALSE}, or
    \item \eqn{{\bf{Y}} = {\bf{K}}^T{\bf{X}}}, if \code{side = "left"} and \code{transposed = TRUE}, or
    \item \eqn{{\bf{Y}} = {\bf{XK}}}, if \code{side = "right"} and \code{transposed = FALSE}, or
    \item \eqn{{\bf{Y}} = {\bf{XK}}^T}, if \code{side = "right"} and \code{transposed = TRUE},
  }
  where \eqn{{\bf{K}}} is the commutation matrix of order \eqn{mn}. The main aim of
  \code{comm.prod} is to do this matrix multiplication \bold{without forming} the
  commutation matrix.
}
\usage{
  comm.prod(m = 1, n = m, x = NULL, transposed = FALSE, side = "left")
}
\arguments{
  \item{m}{ a positive integer row dimension.}
  \item{n}{ a positive integer column dimension.}
  \item{x}{ numeric matrix (or vector).}
  \item{transposed}{ logical. Commutation matrix should be transposed?}
  \item{side}{ a string selecting if commutation matrix is pre-multiplying \code{x}, that is
  \code{side = "left"} or post-multiplying \code{x}, by using \code{side = "right"}.}
}
\details{
  Underlying \code{Fortran} code only uses information provided by \code{\link{comm.info}}
  to performs the matrix multiplication. The commutation matrix is \bold{never} created.
}
\seealso{
  \code{\link{commutation}}
}
\examples{
K42 <- commutation(m = 4, n = 2, matrix = TRUE)
x <- matrix(1:24, ncol = 3)
y <- K42 \%*\% x

z <- comm.prod(m = 4, n = 2, x) # K42 is not stored
all(z == y) # matrices y and z are equal!
}
\keyword{array}
\keyword{algebra}
