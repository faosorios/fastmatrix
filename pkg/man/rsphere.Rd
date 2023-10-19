\name{rsphere}
\alias{rsphere}

\title{Generation of deviates uniformly located on a spherical surface} 
\description{
  Random vector generation uniformly on the sphere.
}
\usage{
rsphere(n = 1, p = 2)
}

\arguments{
  \item{n}{ the number of samples requested}
  \item{p}{ dimension of the unitary sphere}
}

\value{
  If \code{n = 1} a \code{p}-dimensional vector, otherwise a matrix of \code{n}
  rows of random vectors.
}

\details{
  The function \code{rsphere} is an interface to C routines, which make calls to
  subroutines from BLAS. 
}

\references{
  Devroye, L. (1986).
  \emph{Non-Uniform Random Variate Generation}.
  Springer-Verlag, New York.
}

\seealso{
  \code{\link{runif}}
}

\examples{
# generate the sample
z <- rsphere(n = 200)

# scatterplot of a random sample of 200 points uniformly distributed
# on the unit circle
par(pty = "s")
plot(z, xlab = "x", ylab = "y")
title("200 points on the circle", font.main = 1)
}
\keyword{distribution}
\keyword{multivariate}
