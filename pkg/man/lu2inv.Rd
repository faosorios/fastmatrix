\name{lu2inv}
\alias{lu2inv}
\title{Inverse from LU factorization}
\description{
  Invert a square matrix from its LU factorization.
}
\usage{
lu2inv(x)
}
\arguments{
  \item{x}{object representing an LU factorization. This will typically have
  come from a previous call to \code{\link{lu}}.}
}
\value{
  The inverse of the matrix whose LU factorization was given.

  Unsuccessful results from the underlying LAPACK code will result in an
  error giving a positive error code: these can only be interpreted by
  detailed study of the \code{Fortran} code.
}
\source{
  This is an interface to the LAPACK routine \code{DGETRI}. LAPACK is from
  \url{https://netlib.org/lapack/} and its guide is listed in the references.
}
\references{
  Anderson. E., Bai, Z., Bischof, C., Blackford, S., Demmel, J., Dongarra, J.,
  Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A. Sorensen, D. (1999).
  \emph{LAPACK Users' Guide}, 3rd Edition. SIAM.

  Golub, G.H., Van Loan, C.F. (1996).
  \emph{Matrix Computations}, 3rd Edition.
  John Hopkins University Press.
}
\seealso{
  \code{\link{lu}}, \code{\link{solve}}.
}
\examples{
a <- matrix(c(3,2,6,17,4,18,10,-2,-12), ncol = 3)
z <- lu(a)
a \%*\% lu2inv(z)
}
\keyword{algebra}
\keyword{array}
