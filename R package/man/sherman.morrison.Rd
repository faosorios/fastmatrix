\name{sherman.morrison}
\alias{sherman.morrison}
\title{Sherman-Morrison formula}
\description{
  The Sherman-Morrison formula gives a convenient expression for the inverse of the
  rank 1 update \eqn{(\bold{A} + \bold{bd}^T)} where \eqn{\bold{A}} is a \eqn{n\times n}
  matrix and \eqn{\bold{b}}, \eqn{\bold{d}} are \eqn{n}-dimensional vectors. Thus
  \deqn{(\bold{A} + \bold{bd}^T)^{-1} = \bold{A}^{-1} - \frac{\bold{A}^{-1}\bold{bd}^T
  \bold{A}^{-1}}{1 + \bold{d}^T\bold{A}^{-1}\bold{b}}.}
}
\usage{
sherman.morrison(a, b, d = b, inverted = FALSE)
}
\arguments{
  \item{a}{ a numeric matrix.}
  \item{b}{ a numeric vector.}
  \item{d}{ a numeric vector.}
  \item{inverted}{ logical. If \code{TRUE}, \code{a} is supposed to contain its \emph{inverse}.}
}
\details{
  Method of \code{sherman.morrison} calls BLAS level 2 subroutines \code{DGEMV} and
  \code{DGER} for computational efficiency.
}
\value{
  a square matrix of the same order as \code{a}.
}
\examples{
n <- 10
ones <- rep(1, n)
a <- 0.5 * diag(n)
z <- sherman.morrison(a, ones, 0.5 * ones)
z
}
\keyword{array}
\keyword{algebra}
