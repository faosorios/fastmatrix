\name{rank1.update}
\alias{rank1.update}
\title{Rank-one update}
\description{
  Let \eqn{\bold{A}} be a \eqn{n\times n} matrix and \eqn{\bold{u}} an \eqn{n}-dimensional 
  vector. Thus, the transformation
  \deqn{\bold{A} + \alpha\bold{uu}^T,}
  is called a \emph{rank-one update} to \eqn{\bold{A}}, where \eqn{\alpha} is an scalar.
}
\usage{
rank1.update(a, alpha, u, diag = FALSE)
}
\arguments{
  \item{a}{ a numeric matrix.}
  \item{alpha}{ a numeric scalar.}
  \item{u}{ a numeric vector.}
  \item{diag}{ logical. If \code{TRUE}, \code{a} is supposed to be a diagonal matrix and only 
    its diagonal (or a vector) is required.}
}
\details{
  If \code{diag = FALSE}, method of \code{rank1.update} calls BLAS level 2 subroutine \code{DGER} 
  for computational efficiency, otherwise the special structure of the problem is used.
}
\value{
  a square matrix of the same order as \code{a}.
}
\examples{
n <- 10
ones <- rep(1, n)
z <- rank1.update(0.5 * ones, 0.5, ones, diag = TRUE)
z
}
\keyword{array}
\keyword{algebra}
