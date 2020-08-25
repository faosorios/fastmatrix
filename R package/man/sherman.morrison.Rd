\name{sherman.morrison}
\alias{sherman.morrison}
\title{Sherman-Morrison formula}
\description{
  The Sherman-Morrison formula gives a convenient expression for the inverse of the
  rank 1 update \eqn{({\bf{A}} + {\bf{bd}}^T)} where \eqn{{\bf{A}}} is a \eqn{n\times n}
  matrix and \eqn{{\bf{b}}}, \eqn{{\bf{d}}} are \eqn{n}-dimensional vectors. Thus
  \deqn{({\bf{A}} + {\bf{bd}}^T)^{-1} = {\bf{A}}^{-1} - \frac{{\bf{A}}^{-1}{\bf{bd}}^T
  {\bf{A}}^{-1}}{1 + {\bf{d}}^T{\bf{A}}^{-1}{\bf{b}}}.}
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
  Method of \code{sherman.morrison} calls BLAS level 2 subroutines \code{dgemv} and
  \code{dger} for computational efficiency.
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
