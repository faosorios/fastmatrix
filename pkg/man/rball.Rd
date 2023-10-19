\name{rball}
\alias{rball}

\title{Generation of deviates uniformly distributed in a unitary ball} 
\description{
  Random vector generation uniformly in the unitary ball.
}
\usage{
rball(n = 1, p = 2)
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
  The function \code{rball} is an interface to C routines, which make calls to
  subroutines from BLAS. 
}

\references{
  Hormann, W., Leydold, J., Derflinger, G. (2004).
  \emph{Automatic Nonuniform Random Variate Generation}.
  Springer, New York.
}

\seealso{
  \code{\link{runif}}
}

\examples{
# generate the sample
z <- rball(n = 500)

# scatterplot of a random sample of 500 points uniformly distributed
# in the unitary ball
par(pty = "s")
plot(z, xlab = "x", ylab = "y")
title("500 points in the ball", font.main = 1)
}
\keyword{distribution}
\keyword{multivariate}
