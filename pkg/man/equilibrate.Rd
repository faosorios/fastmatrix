\name{equilibrate}
\alias{equilibrate}
\title{Equilibration of a rectangular or symmetric matrix}
\description{
  Equilibrate a rectangular or symmetric matrix using 2-norm.
}
\usage{
equilibrate(x, scale = TRUE)
}
\arguments{
  \item{x}{ a numeric matrix.}
  \item{scale}{ a logical value, \code{x} must be scaled to norm unity?}
}
\value{
  For \code{scale = TRUE}, the equilibrated matrix. The scalings and an approximation
  of the condition number, are returned as attributes \code{"scales"} and \code{"condition"}.
  If \eqn{\bold{X}} is a rectangular matrix, only the columns are equilibrated.
}
\examples{
x <- matrix(c(1, 1, 1,
              1, 2, 1,
              1, 3, 1,
              1, 1,-1,
              1, 2,-1,
              1, 3,-1), ncol = 3, byrow = TRUE)
z <- equilibrate(x)
apply(z, 2, function(x) sum(x^2)) # all 1

xx <- crossprod(x)
equilibrate(xx)
}
\keyword{array}
\keyword{algebra}
