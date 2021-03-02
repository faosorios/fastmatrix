\name{minkowski}
\alias{minkowski}
\title{Computes the p-norm of a vector}
\description{
  Computes a \code{p}-norm of vector \code{x}.  The norm can be the one (\code{p = 1})
  norm, Euclidean (\code{p = 2}) norm, the infinity (\code{p = Inf}) norm. The underlying
  \code{C} or \code{Fortran} code is inspired on ideas of BLAS Level 1.
}
\usage{
minkowski(x, p = 2)
}
\arguments{
  \item{x}{ a numeric vector.}
  \item{p}{ a number, specifying the \emph{type} of norm desired. Possible values include
  real number greater or equal to 1, or Inf, Default value is \code{p = 2}.}
}
\details{
  Method of \code{minkowski} for \code{p = Inf} calls \code{idamax} BLAS function.
  For other values, \code{C} or \code{Fortran} subroutines using unrolled cycles are
  called.
}
\value{
  The vector \code{p}-norm, a non-negative number.
}
\examples{
# a tiny example
x <- rnorm(1000)
minkowski(x, p = 1)
minkowski(x, p = 1.5)
minkowski(x, p = 2)
minkowski(x, p = Inf)

x <- x / minkowski(x)
minkowski(x, p = 2) # equal to 1
}
\keyword{math}
