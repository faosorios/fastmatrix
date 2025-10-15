\name{matrix.sqrt}
\alias{matrix.sqrt}
\title{Matrix square root}
\description{
  This function computes a square root of an \eqn{n\times n} matrix \eqn{\bold{A}}. 
}
\usage{
matrix.sqrt(a, method = "DB", maxiter = 50, tol = 1e-8)
}
\arguments{
  \item{a}{ a square matrix.}
  \item{method}{ the procedure used to obtain the square root. If \code{method = "DB"}
  (the default) the matrix square root is obtained using a Newton's method. 
  If \code{method = "schur"} the Schur decomposition is considered. }
  \item{maxiter}{ the maximum number of iterations. Defaults to \code{50}}
  \item{tol}{ a numeric tolerance.}
}
\details{
  A square root of a square matrix \eqn{\bold{A}} is obtained by solving the 
  equation \eqn{\bold{X}^2 = \bold{A}}, considering the Newton iteration proposed 
  by Denman and Beavers (1976), or alternatively is based on the Schur decomposition.
}
\references{
  Denman, E.D., Beavers, A.N. (1976).
  The matrix sign function and computations in systems.
  \emph{Applied Mathematics and Computation} \bold{2}, 63-94.

  Higham, N.J. (1986).
  Newton's method for the matrix square root.
  \emph{Mathematics of Computation} \bold{46}, 537-549.

  Higham, N.J. (1986).
  \emph{Functions of Matrices: Theory and Computation}.
  Society for Industrial and Applied Mathematics, Philadelphia.
}
\examples{
a <- matrix(c(35,17,3,17,46,11,3,11,12), ncol = 3)
root <- matrix.sqrt(a) # 8 iterations

# just checking
root \%*\% root

root <- matrix.sqrt(a, method = "schur")
}
\keyword{array}
\keyword{algebra}
