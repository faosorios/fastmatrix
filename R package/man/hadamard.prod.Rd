\name{hadamard}
\alias{hadamard}
\title{Hadamard product of two matrices}
\description{
  This function returns the Hadamard or element-wise product of two matrices \code{x}
  and \code{y}, that have the same dimensions.
}
\usage{
  hadamard(x, y = x)
}
\arguments{
  \item{x}{ a numeric matrix or vector.}
  \item{y}{ a numeric matrix or vector.}
}
\value{
  A matrix with the same dimension of \code{x} (and \code{y}) which corresponds to the
  element-by-element product of the two matrices.
}
\references{
  Styan, G.P.H. (1973).
  Hadamard products and multivariate statistical analysis,
  \emph{Linear Algebra and Its Applications} \bold{6}, 217-240.
}
\examples{
x <- matrix(rep(1:10, times = 5), ncol = 5)
y <- matrix(rep(1:5, each = 10), ncol = 5)
z <- hadamard(x, y)
z
}
\keyword{array}
\keyword{algebra}
