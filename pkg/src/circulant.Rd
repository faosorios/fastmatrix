\name{circulant}
\alias{circulant}
\title{Form a symmetric circulant matrix}
\usage{
circulant(x)
}
\arguments{
  \item{x}{the first column to form the circulant matrix.}
}
\description{
  Forms a symmetric circulant matrix using a backwards shift of its first column
}
\value{
  A symmetric circulant matrix.
}
\examples{
x <- c(2,3,5,7,11,13)
circulant(x)
}
\keyword{array}
\keyword{algebra}
