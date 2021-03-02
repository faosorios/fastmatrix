\name{cg}
\alias{cg}
\title{Solve linear systems using the conjugate gradients method}
\usage{
cg(a, b, maxiter = 200, tol = 1e-7)
}
\arguments{
  \item{a}{ a symmetric positive definite matrix containing the coefficients of
  the linear system.}
  \item{b}{ a vector of right-hand sides of the linear system.}
  \item{maxiter}{ the maximum number of iterations. Defaults to \code{200}}
  \item{tol}{ tolerance level for stopping iterations.}
}
\description{
  Conjugate gradients (CG) method is an iterative algorithm for solving linear systems
  with positive definite coefficient matrices.
}
\section{Warning}{
  The underlying \code{C} code does not check for symmetry nor positive definitiveness.
}
\value{
  a vector with the approximate solution, the iterations performed are returned
  as the attribute 'iterations'.
}
\references{
  Golub, G.H., Van Loan, C.F. (1996).
  \emph{Matrix Computations}, 3rd Edition.
  John Hopkins University Press.

  Hestenes, M.R., Stiefel, E. (1952).
  Methods of conjugate gradients for solving linear equations.
  \emph{Journal of Research of the National Bureau of Standards} \bold{49}, 409-436.
}
\seealso{
  \code{\link{jacobi}}, \code{\link{seidel}}, \code{\link{solve}}
}
\examples{
a <- matrix(c(4,3,0,3,4,-1,0,-1,4), ncol = 3)
b <- c(24,30,-24)
z <- cg(a, b)
z # 3 iterations
}
\keyword{algebra}
\keyword{array}
