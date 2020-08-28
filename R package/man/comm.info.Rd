\name{comm.info}
\alias{comm.info}

\title{Compact information to construct the commutation matrix}
\description{
  This function provides the minimum information required to create the commutation matrix.

  The commutation matrix is a square matrix of order \eqn{mn} that, for an \eqn{m\times n}
  matrix \eqn{{\bf{A}}}, transform \code{vec}\eqn{({\bf{A}})} to \code{vec}\eqn{({\bf{A}}^T)}.
}
\usage{
  comm.info(m = 1, n = m, condensed = TRUE)
}
\arguments{
  \item{m}{ a positive integer row dimension.}
  \item{n}{ a positive integer column dimension.}
  \item{condensed}{ logical. Information should be returned in compact form?}
}
\value{
  A list containing the following elements:
  \item{row}{ vector of indexes, each entry represents the row index of the commutation
  matrix.}
  \item{col}{ vector of indexes, each entry represents the column index of the commutation
  matrix. Only present if \code{condensed = FALSE}.}
  \item{m}{ positive integer, row dimension.}
  \item{n}{ positive integer, column dimension.}
}
\details{
  This function returns a list containing two vectors that represent an element of
  the commutation matrix and is accesed by the indexes in vectors \code{row} and \code{col}.
  This information is used by function \code{\link{comm.prod}} to do some operations
  involving the commutation matrix without forming it. This information also can be
  obtained using function \code{\link{commutation}}.
}
\references{
  Magnus, J.R., Neudecker, H. (1979).
  The commutation matrix: some properties and applications.
  \emph{The Annals of Statistics} \bold{7}, 381-394.
}
\seealso{
  \code{\link{commutation}}, \code{\link{comm.prod}}
}
\examples{
z <- comm.info(m = 3, n = 2, condensed = FALSE)
z # where are the ones in commutation matrix of order '3,2'?

K32 <- commutation(m = 3, n = 2, matrix = TRUE)
K32 # only recommended if m and n are very small
}
\keyword{array}
