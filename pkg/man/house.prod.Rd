\name{house.prod}
\alias{house.prod}
\title{Applying Householder matrices}
\description{
  This function applies the Householder matrix defined by a nonzero vector \eqn{\bold{x}} to 
  a matrix \eqn{\bold{A}}. Thus, the output of this function is the matrix-matrix operation:
  \itemize{
    \item \eqn{\bold{Y} = \bold{PA}}, if \code{side = "left"}, or
    \item \eqn{\bold{Y} = \bold{AP}}, if \code{side = "right"}.
  }
}
\usage{
house.prod(a, x, side = "left")
}
\arguments{
  \item{a}{ a numeric rectangular matrix.}
  \item{x}{ an \eqn{n}-dimensional numeric vector.}
  \item{side}{ a string selecting if Householder matrix is pre-multiplying \code{a}, that is
  \code{side = "left"} or post-multiplying \code{a}, by using \code{side = "right"}.}
}
\details{
  Underlying code \bold{never} entail the explicit formation of the Householder matrix.
}
\references{
  Golub, G.H., van Loan, C.F. (1996).
  \emph{Matrix Computations}, 3rd Edition. 
  The Johns Hopkins University Press, Baltimore.
}
\examples{
x <- c(3,1,5,1)
u <- house(x)$u
house.prod(u, x)

y <- as.matrix(x)
house.prod(y, x)
}
\keyword{array}
\keyword{algebra}
