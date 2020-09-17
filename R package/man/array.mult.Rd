\name{array.mult}
\alias{array.mult}
\title{Array multiplication}
\description{
  Multiplication of 3-dimensional arrays was first introduced by Bates and Watts (1980).
  More extensions and technical details can be found in Wei (1998).
}
\usage{
array.mult(a, b, x)
}
\arguments{
  \item{a}{ a numeric matrix.}
  \item{b}{ a numeric matrix.}
  \item{x}{ a three-dimensional array.}
}
\details{
  Let \eqn{\bold{X} = (x_{tij})} be a 3-dimensional \eqn{n\times p\times q} where
  indices \eqn{t, i} and \eqn{j} indicate face, row and column, respectively. The
  product \eqn{\bold{Y} = \bold{AXB}} is an \eqn{n\times r\times s} array, with
  \eqn{\bold{A}} and \eqn{\bold{B}} are \eqn{r\times p} and \eqn{q\times s} matrices
  respectively. The elements of \eqn{\bold{Y}} are defined as:
  \deqn{y_{tkl} = \sum\limits_{i=1}^p\sum\limits_{j=1}^q a_{ki}x_{tij}b_{jl}}
}
\value{
  \code{array.mult} returns a 3-dimensional array of dimension \eqn{n\times r\times s}.
}
\references{
  Bates, D.M., Watts, D.G. (1980).
  Relative curvature measures of nonlinearity.
  \emph{Journal of the Royal Statistical Society, Series B} \bold{42}, 1-25.

  Wei, B.C. (1998).
  \emph{Exponential Family Nonlinear Models}.
  Springer, New York.
}
\seealso{
  \code{\link{array}}, \code{\link{matrix}}, \code{\link{bracket.prod}}.
}
\examples{
x <- array(0, dim = c(2,3,3)) # 2 x 3 x 3 array
x[,,1] <- c(1,2,2,4,3,6)
x[,,2] <- c(2,4,4,8,6,12)
x[,,3] <- c(3,6,6,12,9,18)

a <- matrix(1, nrow = 2, ncol = 3)
b <- matrix(1, nrow = 3, ncol = 2)

y <- array.mult(a, b, x) # a 2 x 2 x 2 array
y
}
\keyword{array}
\keyword{algebra}
