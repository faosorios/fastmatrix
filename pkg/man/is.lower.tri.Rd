\name{is.lower.tri}
\alias{is.lower.tri}
\alias{is.upper.tri}
\title{Check if a matrix is lower or upper triangular}
\description{
  Returns \code{TRUE} if the given matrix is lower or upper triangular matrix.
}
\usage{
is.lower.tri(x, diag = FALSE)
is.upper.tri(x, diag = FALSE)
}
\arguments{
  \item{x}{a matrix of other \R object with \code{length(dim(x)) == 2}.}
  \item{diag}{logical. Should the diagonal be included?}
}
\value{
  Check if a matrix is lower or upper triangular. You can also include diagonal to the check.
}
\seealso{
  \code{\link{lower.tri}}, \code{\link{upper.tri}}
}
\examples{
  x <- matrix(rnorm(10 * 3), ncol = 3)
  R <- chol(crossprod(x))

  is.lower.tri(R)
  is.upper.tri(R)
}
\keyword{array}
