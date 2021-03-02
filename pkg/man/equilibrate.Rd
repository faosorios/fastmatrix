\name{equilibrate}
\alias{equilibrate}
\title{Column equilibration of a rectangular matrix}
\description{
  Equilibrate the columns of a rectangular matrix using 2-norm.
}
\usage{
equilibrate(x, scale = TRUE)
}
\arguments{
  \item{x}{ a numeric matrix.}
  \item{scale}{ a logical value, the columns of \code{x} must be scaled to norm unity?}
}
\value{
  For \code{scale = TRUE}, the equilibrated (each column scaled to norm one) matrix.
  The scalings and an approximation of the reciprocal condition number, are returned
  as attributes \code{"scales"} and \code{"condition"}.
}
\examples{
x <- matrix(c(1, 1, 1,
              1, 2, 1,
              1, 3, 1,
              1, 1,-1,
              1, 2,-1,
              1, 3,-1), ncol = 3, byrow = TRUE)
x <- equilibrate(x)
apply(x, 2, function(x) sum(x^2)) # all 1
}
\keyword{array}
\keyword{algebra}
