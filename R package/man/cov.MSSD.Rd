\name{cov.MSSD}
\alias{cov.MSSD}
\title{Variance and covariance matrices}
\description{
  Returns a list containing the mean and covariance matrix of the data.
}
\usage{
cov.MSSD(x)
}
\arguments{
  \item{x}{a matrix or data frame. As usual, rows are observations and columns are
    variables.}
}
\value{
  A list containing the following named components:
  \item{mean}{an estimate for the center (mean) of the data.}
  \item{cov}{the estimated covariance matrix.}
}
\details{
  This procedure uses the Holmes-Mergen method using the difference between each
  successive pairs of observations also known as Mean Square Successive Method (MSSD)
  to estimate the covariance matrix.
}
\seealso{
  \code{\link{cov}} and \code{\link{var}}.
}
\references{
  Holmes, D.S., Mergen, A.E. (1993).
  Improving the performance of the \eqn{T^2} control chart.
  \emph{Quality Engineering} \bold{5}, 619-625.
}
\examples{
xy <- cbind(x = 1:10, y = c(1:3, 8:5, 8:10))
z0 <- cov(xy)
z0
z1 <- cov.MSSD(xy)
z1
}
\keyword{multivariate}
