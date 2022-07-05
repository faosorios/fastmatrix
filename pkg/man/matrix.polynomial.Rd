\name{matrix.polynomial}
\alias{matrix.polynomial}
\title{Evaluates a real general matrix polynomial}
\usage{
matrix.polynomial(a, coef = rep(1, power + 1), power = length(coef))
}
\arguments{
  \item{a}{ a numeric square matrix of order \eqn{n} by \eqn{n} for which the polinomial
  is to be computed.}
  \item{coef}{ numeric vector containing the coefficients of the polinomial in order of
  increasing power.}
  \item{power}{ a numeric exponent (which is forced to be an integer). If the exponent
  is zero, a multiple of the identity matrix is returned.}
}
\description{
  Given \eqn{c_0,c_1,\dots,c_n} coefficients of the polynomial and \eqn{\bold{A}}
  an \eqn{n} by \eqn{n} matrix. This function computes the matrix polynomial
  \deqn{\bold{B} = \sum\limits_{k=0}^n c_k\bold{A}^k,}
  using Horner's scheme, where \eqn{\bold{A}^0 = \bold{I}} is the \eqn{n} by \eqn{n} identity matrix.
}
\value{
  Returns an \eqn{n} by \eqn{n} matrix.
}

\examples{
a <- matrix(c(1, 3, 2, -5, 1, 7, 1, 5, -4), ncol = 3, byrow = TRUE)
cf <- c(3, 1, 2)
b <- matrix.polynomial(a, cf)
b # 3 * diag(3) + a + 2 * a %*% a
b <- matrix.polynomial(a, power = 2)
b # diag(3) + a + a %*% a
}
\keyword{array}
