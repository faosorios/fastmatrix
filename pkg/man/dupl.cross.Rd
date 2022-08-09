\name{dupl.cross}
\alias{dupl.cross}

\title{Matrix crossproduct envolving the duplication matrix}
\description{
  Given the order of two duplication matrices and a conformable matrix \eqn{\bold{X}},
  this function performs the operation: \eqn{\bold{Y} = \bold{D}_n^T\bold{X}\bold{D}_k},
  where \eqn{\bold{D}_n} and \eqn{\bold{D}_k} are duplication matrices of order \eqn{n}
  and \eqn{k}, respectively.
}
\usage{
dupl.cross(n = 1, k = n, x = NULL)
}
\arguments{
  \item{n}{ order of the duplication matrix used pre-multiplying \code{x}.}
  \item{k}{ order of the duplication matrix used post-multiplying \code{x}. By default
  \code{k = n} is used.}
  \item{x}{ numeric matrix, this argument is required.}
}
\details{
  This function calls \code{\link{dupl.prod}} to performs the matrix multiplications required
  but \bold{without forming} any duplication matrices.
}
\seealso{
  \code{\link{dupl.prod}}
}
\examples{
D2 <- duplication(n = 2, matrix = TRUE)
D3 <- duplication(n = 3, matrix = TRUE)
x <- matrix(1, nrow = 9, ncol = 4)
y <- t(D3) \%*\% x \%*\% D2

z <- dupl.cross(n = 3, k = 2, x) # D2 and D3 are not stored
all(z == y) # matrices y and z are equal!

x <- matrix(1, nrow = 9, ncol = 9)
z <- dupl.cross(n = 3, x = x) # same matrix is used to pre- and post-multiplying x
z # print result
}
\keyword{array}
\keyword{algebra}
