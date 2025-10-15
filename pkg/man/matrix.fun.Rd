\name{matrix.fun}
\alias{matrix.fun}
\title{Evaluate a matrix function}
\description{
  This function computes the matrix function \eqn{\bold{F} =f(\bold{A})} where 
  \eqn{\bold{A}} is upper triangular by applying a Parlett recurrence.
}
\usage{
matrix.fun(a, FUN = "log")
}
\arguments{
  \item{a}{ an upper triangular matrix.}
  \item{FUN}{ the function to be applied, by default \code{"log"}.}
}
\details{
  The used-defined function \code{FUN} is evaluated at the triangular matrix 
  argument. This function can be used in conjunction with Schur decomposition 
  to evaluate the function of a matrix.
}
\references{
  Higham, N.J. (1986).
  \emph{Functions of Matrices: Theory and Computation}.
  Society for Industrial and Applied Mathematics, Philadelphia.
}
\examples{
a <- matrix(c(1,2,3,0,3,4,0,0,5), ncol = 3, byrow = TRUE)
fnc <- function(x) (1 + x) / x
f <- matrix.fun(a, FUN = fnc)
f

a <- matrix(c(-49,24,-64,31), ncol = 2, byrow = TRUE)
z <- schur(a)
m <- z$m
u <- z$vectors
m <- matrix.fun(m, FUN = exp)
u \%*\% m \%*\% t(u) # exp(a)
}
\keyword{array}
\keyword{algebra}
