\name{power.method}
\alias{power.method}
\title{Power method to approximate dominant eigenvalue and eigenvector}
\description{
  The power method seeks to determine the eigenvalue of maximum modulus, and a corresponding
  eigenvector.
}
\usage{
  power.method(x, only.value = FALSE, maxiter = 100, tol = 1e-8)
}
\arguments{
  \item{x}{ a symmetric matrix.}
  \item{only.value}{ if \code{TRUE}, only the dominant eigenvalue is returned, otherwise
  both dominant eigenvalue and eigenvector are returned.}
  \item{maxiter}{ the maximum number of iterations. Defaults to \code{100}}
  \item{tol}{ a numeric tolerance.}
}
\value{
  When \code{only.value} is not true, as by default, the result is a list with components
  \code{"value"} and \code{"vector"}. Otherwise only the dominan eigenvalue is returned.
  The performed number of iterations to reach convergence is returned as attribute \code{"iterations"}.
}
\seealso{
  \code{\link{eigen}} for eigenvalues and eigenvectors computation.
}
\examples{
n <- 1000
x <- .5 * diag(n) + 0.5 * matrix(1, nrow = n, ncol = n)

# dominant eigenvalue must be (n + 1) / 2
z <- power.method(x, only.value = TRUE)
}
\keyword{array}
\keyword{algebra}
