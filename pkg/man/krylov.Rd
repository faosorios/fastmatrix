\name{krylov}
\alias{krylov}
\title{Computes a Krylov matrix}
\usage{
krylov(a, b, m = ncol(a))
}
\arguments{
  \item{a}{ a numeric square matrix of order \eqn{n} by \eqn{n} for which the Krylov
  matrix is to be computed.}
  \item{b}{ a numeric vector of length \eqn{n}.}
  \item{m}{ length of the Krylov sequence.}
}
\description{
  Given \eqn{\bold{A}} an \eqn{n} by \eqn{n} real matrix and an \eqn{n}-vector \eqn{\bold{b}},
  this function constructs the Krylov matrix \eqn{\bold{K}}, where
  \deqn{\bold{K} = [\bold{b},\bold{Ab},\dots,\bold{A}^{m-1}\bold{b}].}
}
\value{
  Returns an \eqn{n} by \eqn{m} matrix.
}
\examples{
a <- matrix(c(1, 3, 2, -5, 1, 7, 1, 5, -4), ncol = 3, byrow = TRUE)
b <- c(1, 1, 1)
k <- krylov(a, b, m = 4)
k
}
\keyword{array}
