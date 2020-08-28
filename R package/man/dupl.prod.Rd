\name{dupl.prod}
\alias{dupl.prod}
\title{Matrix multiplication envolving the duplication matrix}
\description{
  Given the order of a duplication and matrix \code{x}, performs one of the matrix-matrix
  operations:
  \itemize{
    \item \eqn{{\bf{Y}} = {\bf{DX}}}, if \code{side = "left"} and \code{transposed = FALSE}, or
    \item \eqn{{\bf{Y}} = {\bf{D}}^T{\bf{X}}}, if \code{side = "left"} and \code{transposed = TRUE}, or
    \item \eqn{{\bf{Y}} = {\bf{XD}}}, if \code{side = "right"} and \code{transposed = FALSE}, or
    \item \eqn{{\bf{Y}} = {\bf{XD}}^T}, if \code{side = "right"} and \code{transposed = TRUE},
  }
  where \eqn{{\bf{D}}} is the duplication matrix of order \eqn{n}. The main aim of
  \code{dupl.prod} is to do this matrix multiplication \bold{without forming} the
  duplication matrix.
}
\usage{
  dupl.prod(n = 1, x, transposed = FALSE, side = "left")
}
\arguments{
  \item{n}{ order of the duplication matrix.}
  \item{x}{ numeric matrix (or vector).}
  \item{transposed}{ logical. Duplication matrix should be transposed?}
  \item{side}{ a string selecting if duplication matrix is pre-multiplying \code{x}, that is
  \code{side = "left"} or post-multiplying \code{x}, by using \code{side = "right"}.}
}
\details{
  Underlying \code{C} code only uses information provided by \code{\link{dupl.info}} to
  performs the matrix multiplication. The duplication matrix is \bold{never} created.
}
\seealso{
  \code{\link{duplication}}
}
\examples{
D4 <- duplication(n = 4, matrix = TRUE)
x <- matrix(1, nrow = 16, ncol = 2)
y <- crossprod(D4, x)

z <- dupl.prod(n = 4, x, transposed = TRUE) # D4 is not stored
all(z == y) # matrices y and z are equal!
}
\keyword{array}
\keyword{algebra}
