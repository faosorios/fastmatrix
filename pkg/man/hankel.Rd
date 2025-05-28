\name{hankel}
\alias{hankel}
\title{Form a symmetric Hankel matrix}
\usage{
hankel(x, y = NULL)
}
\arguments{
  \item{x}{the first column to form the Hankel matrix.}
  \item{y}{the last column of the Hankel matrix. If \code{y} is not provided only 
  its first \eqn{n} skew diagonals are formed and the remaining elements are zeros. 
  Otherwise, it is assumed that \eqn{x_n = y_1}, and the first entry of vector \code{y} 
  is discarded.}
}
\description{
  Forms a symmetric Hankel matrix of order \eqn{n} from the values in vector \eqn{\bold{x}} 
  and optionally the vector \eqn{\bold{y}}.
}
\value{
  A symmetric Hankel matrix of order \eqn{n}.
}
\examples{
x <- 1:4
y <- c(4,6,8,10)

# H4
hankel(x)

# H({1,2,3,4},{4,6,8,10})
hankel(x, y)
}
\keyword{array}
\keyword{algebra}
