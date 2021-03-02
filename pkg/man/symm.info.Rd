\name{symm.info}
\alias{symm.info}

\title{Compact information to construct the symmetrizer matrix}
\description{
  This function provides the information required to create the symmetrizer matrix.
}
\usage{
symm.info(n = 1)
}
\arguments{
  \item{n}{ order of the symmetrizer matrix.}
}
\value{
  A list containing the following elements:
  \item{row}{ vector of indexes, each entry represents the row index of the symmetrizer
  matrix.}
  \item{col}{ vector of indexes, each entry represents the column index of the symmetrizer
  matrix.}
  \item{val}{ vector of values, each entry represents the value of the symmetrizer matrix
  at element given by \code{row} and \code{col} indexes.}
  \item{order}{ order of the symmetrizer matrix.}
}
\details{
  This function returns a list containing vectors that represent an element of the
  symmetrizer matrix and is accesed by the indexes in vectors \code{row}, \code{col}
  and values contained in \code{val}. This information is used by function \code{\link{symm.prod}}
  to do some operations involving the symmetrizer matrix without forming it. This
  information also can be obtained using function \code{\link{symmetrizer}}.
}
\seealso{
  \code{\link{symmetrizer}}, \code{\link{symm.prod}}
}
\examples{
z <- symm.info(n = 3)
z # elements in symmetrizer matrix of order 3

N3 <- symmetrizer(n = 3, matrix = TRUE)
N3 # only recommended if n is very small
}
\keyword{array}
