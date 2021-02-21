\name{seidel}
\alias{seidel}
\title{Solve linear systems using the Gauss-Seidel method}
\usage{
seidel(a, b, start, maxiter = 200, tol = 1e-7)
}
\arguments{
  \item{a}{ a square numeric matrix containing the coefficients of the linear system.}
  \item{b}{ a vector of right-hand sides of the linear system.}
  \item{start}{ a vector for initial starting point.}
  \item{maxiter}{ the maximum number of iterations. Defaults to \code{200}}
  \item{tol}{ tolerance level for stopping iterations.}
}
\description{
  Gauss-Seidel method is an iterative algorithm for solving a system of linear equations.
}
\details{
  Let \eqn{\bold{D}}, \eqn{\bold{L}}, and \eqn{\bold{U}} denote the diagonal, lower
  triangular and upper triangular parts of a matrix \eqn{\bold{A}}. Gauss-Seidel method
  solve the equation \eqn{\bold{Ax} = \bold{b}}, iteratively by rewriting \eqn{(\bold{L}
  + \bold{D})\bold{x} + \bold{Ux} = \bold{b}}. Assuming that \eqn{\bold{L} + \bold{D}} is
  nonsingular leads to the iteration formula
  \deqn{\bold{x}^{(k+1)} = -(\bold{L} + \bold{D})^{-1}\bold{U}\bold{x}^{(k)} + (\bold{L}
  + \bold{D})^{-1}\bold{b}}
}
\value{
  a vector with the approximate solution, the iterations performed are returned
  as the attribute 'iterations'.
}
\references{
  Golub, G.H., Van Loan, C.F. (1996).
  \emph{Matrix Computations}, 3rd Edition.
  John Hopkins University Press.
}
\seealso{
  \code{\link{jacobi}}
}
\examples{
a <- matrix(c(5,-3,2,-2,9,-1,3,1,-7), ncol = 3)
b <- c(-1,2,3)
start <- c(1,1,1)
z <- seidel(a, b, start)
z # 10 iterations
}
\keyword{algebra}
\keyword{array}
