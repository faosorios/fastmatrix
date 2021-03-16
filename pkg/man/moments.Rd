\name{moments}
\alias{moments}
\title{Central moments}
\description{
  It calculates up to fourth central moments (or moments about the mean), and the
  skewness and kurtosis coefficients using an online algorithm.
}
\usage{
moments(x)
}
\arguments{
  \item{x}{a numeric vector containing the sample observations.}
}
\value{
  A list containing \code{second}, \code{third} and \code{fourth} central moments,
  and \code{skewness} and \code{kurtosis} coefficients.
}
\details{
  The \eqn{k}-th central moment is defined as
  \deqn{m_k = \frac{1}{n}\sum_{i=1}^n (x_i - \overline{x})^k.}
  In particular, the second central moment is the variance of the sample. The sample
  skewness and kurtosis are defined, respectively, as
  \deqn{b_1 = \frac{m_3}{s^3}, \qquad b_2 = \frac{m_4}{s^4} - 3,}
  where \eqn{s} denotes de standard deviation.
}
\seealso{
  \code{\link{var}}.
}
\references{
  Spicer, C.C. (1972).
  Algorithm AS 52: Calculation of power sums of deviations about the mean.
  \emph{Applied Statistics} \bold{21}, 226-227.
}
\examples{
set.seed(149)
x <- rnorm(1000)
z <- moments(x)
z
}
\keyword{univar}
