\name{scaled.condition}
\alias{scaled.condition}
\title{Scaled condition number}
\description{
  Compute the scaled condition number of a rectangular matrix.
}
\usage{
scaled.condition(x, scales = FALSE)
}
\arguments{
  \item{x}{ a numeric rectangular matrix.}
  \item{scales}{ a logical value indicating whether the scaling factors that allow balancing 
  the columns of \code{x} should be returned by the function.}
}
\value{
  The columns of a rectangular matrix \code{x} are equilibrated (but not centered), 
  then the scaled condition number is computed following the guidelines of Belsley (1990).
  If requested, the column scalings are returned as the attribute \code{'scales'}.
}
\references{
  Belsley, D.A. (1990).
  \emph{Conditioning Diagnostics: Collinearity and Weak Data in Regression}.
  Wiley, New York.
}
\examples{
x <- matrix(c(1, 1, 1,
              1, 2, 1,
              1, 3, 1,
              1, 1,-1,
              1, 2,-1,
              1, 3,-1), ncol = 3, byrow = TRUE)
scaled.condition(x)
}
\keyword{array}
\keyword{algebra}
