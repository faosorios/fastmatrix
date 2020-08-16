\name{sweep.operator}
\alias{sweep.operator}
\title{Gauss-Jordan sweep operator for symmetric matrices}
\description{
  Perform the sweep operation (or reverse sweep) on the \code{k}-th element of the
  diagonal of a symmetric matrix.
}
\usage{
  sweep.operator(x, k = 1, reverse = FALSE)
}
\arguments{
  \item{x}{ a symmetric matrix.}
  \item{k}{ element of the diagonal which will be sweeped.}
  \item{reverse}{ logical. If \code{reverse = TRUE} the reverse sweep is performed.}
}
\details{
  The symmetric sweep operator is a powerful tool in computational statistics with
  uses in stepwise regression, conditional multivariate normal distributions, MANOVA,
  and more.
}
\value{
  a square matrix of the same order as \code{x}.
}
\references{
  Goodnight, J.H. (1979).
  A tutorial on the SWEEP operator.
  \emph{The American Statistician} \bold{33}, 149-158.
}
\examples{
x <- matrix(rnorm(1000 * 100), ncol = 100)
xx <- crossprod(x)
y <- sweep.operator(xx, k = 1)
}
\keyword{array}
\keyword{algebra}
