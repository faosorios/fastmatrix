\name{rmnorm}
\alias{rmnorm}
\title{Multivariate normal random deviates}
\description{
  Random number generation from the multivariate normal (Gaussian) distribution.
}
\usage{
rmnorm(n = 1, mean = rep(0, nrow(Sigma)), Sigma = diag(length(mean)))
}

\arguments{
  \item{n}{ the number of samples requested}
  \item{mean}{ a vector giving the means of each variable}
  \item{Sigma}{ a positive-definite covariance matrix}
}

\value{
  If \eqn{n = 1} a vector of the same length as \code{mean}, otherwise a
  matrix of \eqn{n} rows of random vectors.
}

\details{
  The function \code{rmnorm} is an interface to \code{C} routines, which make calls
  to subroutines from LAPACK. The matrix decomposition is internally done using the
  Cholesky decomposition. If \code{Sigma} is not non-negative definite then there
  will be a warning message.
}

\references{
  Devroye, L. (1986).
  \emph{Non-Uniform Random Variate Generation}.
  Springer-Verlag, New York.
}

\seealso{
  \code{\link{rnorm}}
}

\examples{
# covariance parameters
Sigma <- matrix(c(10,3,3,2), ncol = 2)
Sigma

# generate the sample
y <- rmnorm(n = 1000, Sigma = Sigma)
var(y)

# scatterplot of a random bivariate normal sample with mean
# vector zero and covariance matrix 'Sigma'
par(pty = "s")
plot(y, xlab = "", ylab = "")
title("bivariate normal sample", font.main = 1)

# QQ-plot of transformed distances
z <- wilson.hilferty(y)
par(pty = "s")
qqnorm(z, main = "Transformed distances Q-Q plot")
abline(c(0,1), col = "red", lwd = 2, lty = 2)
}
\keyword{distribution}
\keyword{multivariate}
