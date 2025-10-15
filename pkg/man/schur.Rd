\name{schur}
\alias{schur}
\title{The Schur decomposition of a square matrix}
\usage{
schur(x, vectors = TRUE)
}
\arguments{
  \item{x}{a square numeric matrix to be decomposed.}
  \item{vectors}{ logical, if \code{TRUE} (the default), then Schur vectors are returned.}
}
\description{
  \code{schur} computes the Schur decomposition of an \eqn{n \times n}{n-by-n} real matrix 
  \eqn{\bold{A}}.
}
\details{
  For an \eqn{n \times n}{n-by-n} real matrix \eqn{\bold{A}}, the Schur decomposition 
  is given by,
  \deqn{\bold{A} = \bold{VMV}^T}
  where \eqn{\bold{V}} is an orthogonal matrix and \eqn{\bold{M}} is an upper triangular 
  matrix. The column vectors \eqn{\bold{V}} (if requested) are the Schur vectors of 
  \eqn{\bold{A}}, and \eqn{\bold{M}} is the Schur form of \eqn{\bold{A}}.

  Unsuccessful results from the underlying LAPACK code will result in an error giving 
  a error code: these can only be interpreted by detailed study of the \code{Fortran} code.
}
\value{
  The Schur decomposition of the matrix as computed by LAPACK. The components in
  the returned value correspond directly to the values returned by \code{DGEES}.
  \item{m}{ a matrix with the same dimensions as \code{x}. The upper triangle contains 
    the \eqn{\bold{M}} matrix of the decomposition.}
  \item{values}{ a vector containing the \eqn{n} eigenvalues of \code{x}, these values 
    \bold{are not ordered}.}
  \item{vectors}{ an \eqn{n \times n}{n-by-n} matrix whose columns contain the eigenvectors 
    of \code{x}, only available if it is requested.}
}
\references{
  Anderson. E., Bai, Z., Bischof, C., Blackford, S., Demmel, J., Dongarra, J.,
  Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A. Sorensen, D. (1999).
  \emph{LAPACK Users' Guide}, 3rd Edition. SIAM.

  Golub, G.H., Van Loan, C.F. (1996).
  \emph{Matrix Computations}, 3rd Edition.
  John Hopkins University Press.
}
\examples{
a <- matrix(c(7,12,-2,-3), ncol = 2)
z <- schur(a)
z # information of Schur decomposition

x <- matrix(c(0,0,1,2,1,0,2,2,1), ncol = 3)
z <- schur(x)
z # complex eigenvalues
}
\keyword{algebra}
\keyword{array}
