\name{bracket.prod}
\alias{bracket.prod}
\title{Bracket product}
\description{
  Bracket product of a matrix and a 3-dimensional array.
}
\usage{
bracket.prod(a, x)
}
\arguments{
  \item{a}{ a numeric matrix.}
  \item{x}{ a three-dimensional array.}
}
\details{
  Let \eqn{\bold{X} = (x_{tij})} be a 3-dimensional \eqn{n\times p\times q} array and
  \eqn{\bold{A}} an \eqn{m\times n} matrix, then \eqn{\bold{Y} = [\bold{A}][\bold{X}]}
  is called the bracket product of \eqn{\bold{A}} and \eqn{\bold{X}}, that is an \eqn{m
  \times p\times q} with elements
  \deqn{y_{tij} = \sum\limits_{k=1}^n a_{tk}x_{kij}}
}
\value{
  \code{bracket.prod} returns a 3-dimensional array of dimension \eqn{m\times p\times q}.
}
\references{
  Wei, B.C. (1998).
  \emph{Exponential Family Nonlinear Models}.
  Springer, New York.
}
\seealso{
  \code{\link{array}}, \code{\link{matrix}}, \code{\link{array.mult}}.
}
\examples{
x <- array(0, dim = c(2,3,3)) # 2 x 3 x 3 array
x[,,1] <- c(1,2,2,4,3,6)
x[,,2] <- c(2,4,4,8,6,12)
x[,,3] <- c(3,6,6,12,9,18)

a <- matrix(1, nrow = 3, ncol = 2)

y <- bracket.prod(a, x) # a 3 x 3 x 3 array
y
}
\keyword{array}
\keyword{algebra}
