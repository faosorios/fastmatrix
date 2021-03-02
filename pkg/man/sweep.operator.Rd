\name{sweep.operator}
\alias{sweep.operator}
\title{Gauss-Jordan sweep operator for symmetric matrices}
\description{
  Perform the sweep operation (or reverse sweep) on the diagonal elements of a
  symmetric matrix.
}
\usage{
sweep.operator(x, k = 1, reverse = FALSE)
}
\arguments{
  \item{x}{ a symmetric matrix.}
  \item{k}{ elements (if \code{k} is vector) of the diagonal which will be sweeped.}
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
# tiny example of regression, last column contains 'y'
xy <- matrix(c(1, 1, 1, 1,
               1, 2, 1, 3,
               1, 3, 1, 3,
               1, 1,-1, 2,
               1, 2,-1, 2,
               1, 3,-1, 1), ncol = 4, byrow = TRUE)
z <- crossprod(xy)
z <- sweep.operator(z, k = 1:3)
cf <- z[1:3,4] # regression coefficients
RSS <- z[4,4]  # residual sum of squares

# an example not that small
x <- matrix(rnorm(1000 * 100), ncol = 100)
xx <- crossprod(x)
z <- sweep.operator(xx, k = 1)
}
\keyword{array}
\keyword{algebra}
