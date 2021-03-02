\name{lu}
\alias{lu}
\alias{lu.default}
\alias{solve.lu}
\alias{is.lu}
\title{The LU factorization of a square matrix}
\usage{
lu(x)
\method{lu}{default}(x)

\method{solve}{lu}(a, b, \dots)

is.lu(x)
}
\arguments{
  \item{x}{a square numeric matrix whose LU factorization is to be computed.}
  \item{a}{an LU factorization of a square matrix.}
  \item{b}{a vector or matrix of right-hand sides of equations.}
  \item{\dots}{further arguments passed to or from other methods}
}
\description{
  \code{lu} computes the LU factorization of a matrix.
}
\details{
  The LU factorization plays an important role in many numerical procedures. In
  particular it is the basic method to solve the equation \eqn{\bold{Ax} = \bold{b}}
  for given matrix \eqn{\bold{A}}, and vector \eqn{\bold{b}}.

  \code{solve.lu} is the method for \code{\link{solve}} for \code{lu} objects.

  \code{is.lu} returns \code{TRUE} if \code{x} is a \code{\link{list}}
  and \code{\link{inherits}} from \code{"lu"}.% not more checks, for speed.

  Unsuccessful results from the underlying LAPACK code will result in an
  error giving a positive error code: these can only be interpreted by
  detailed study of the \code{Fortran} code.
}
\value{
  The LU factorization of the matrix as computed by LAPACK. The components in
  the returned value correspond directly to the values returned by \code{DGETRF}.
  \item{lu}{a matrix with the same dimensions as \code{x}.
    The upper triangle contains the \eqn{\bold{U}} of the decomposition
    and the strict lower triangle contains information on the \eqn{\bold{L}} of
    the factorization.}
  \item{pivot}{information on the pivoting strategy used during
    the factorization.}
}
\note{
  To compute the determinant of a matrix (do you \emph{really} need it?),
  the LU factorization is much more efficient than using eigenvalues
  (\code{\link{eigen}}).  See \code{\link{det}}.

  LAPACK uses column pivoting and does not attempt to detect rank-deficient matrices.
}
\references{
  Anderson. E., Bai, Z., Bischof, C., Blackford, S., Demmel, J., Dongarra, J.,
  Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A. Sorensen, D. (1999).
  \emph{LAPACK Users' Guide}, 3rd Edition. SIAM.
  (Available at \url{http://www.netlib.org/lapack/lug/lapack_lug.html}).

  Golub, G.H., Van Loan, C.F. (1996).
  \emph{Matrix Computations}, 3rd Edition.
  John Hopkins University Press.
}
\seealso{
  \code{\link{extractL}},  \code{\link{extractU}},  \code{\link{constructX}} for
  reconstruction of the matrices,
  \code{\link{lu2inv}}
}
\examples{
a <- matrix(c(3,2,6,17,4,18,10,-2,-12), ncol = 3)
z <- lu(a)
z # information of LU factorization

# computing det(a)
prod(diag(z$lu)) # product of diagonal elements of U

# solve linear equations
b <- matrix(1:6, ncol = 2)
solve(z, b)
}
\keyword{algebra}
\keyword{array}
