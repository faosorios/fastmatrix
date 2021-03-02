\name{dupl.info}
\alias{dupl.info}

\title{Compact information to construct the duplication matrix}
\description{
  This function provides the minimum information required to create the duplication matrix.
}
\usage{
dupl.info(n = 1, condensed = TRUE)
}
\arguments{
  \item{n}{ order of the duplication matrix.}
  \item{condensed}{ logical. Information should be returned in compact form?}
}
\value{
  A list containing the following elements:
  \item{row}{ vector of indexes, each entry represents the row index of the duplication
  matrix. Only present if \code{condensed = FALSE}.}
  \item{col}{ vector of indexes, each entry represents the column index of the duplication
  matrix.}
  \item{order}{ order of the duplication matrix.}
}
\details{
  This function returns a list containing two vectors that represent an element of
  the duplication matrix and is accesed by the indexes in vectors \code{row} and \code{col}.
  This information is used by function \code{\link{dupl.prod}} to do some operations
  involving the duplication matrix without forming it. This information also can be
  obtained using function \code{\link{duplication}}
}
\seealso{
  \code{\link{duplication}}, \code{\link{dupl.prod}}
}
\examples{
z <- dupl.info(n = 3, condensed = FALSE)
z # where are the ones in duplication of order 3?

D3 <- duplication(n = 3, matrix = TRUE)
D3 # only recommended if n is very small
}
\keyword{array}
