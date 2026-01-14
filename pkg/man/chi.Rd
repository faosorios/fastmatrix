\name{chi}
\alias{chi}
\alias{dchi}
\alias{pchi}
\alias{qchi}
\alias{rchi}
\title{The chi distribution}
\description{
  Density, distribution function, quantile function and random generation for the
  chi distribution with \code{df} degrees of freedom.
}
\usage{
dchi(x, df = 1, log = FALSE)
pchi(q, df = 1, lower.tail = TRUE, log.p = FALSE)
qchi(p, df = 1, lower.tail = TRUE, log.p = FALSE)
rchi(n, df = 1)
}
\arguments{
  \item{x, q}{ vector of quantiles.}
  \item{p}{ vector of probabilities.}
  \item{n}{ number of observations.}
  \item{df}{ degrees of freedom (non-negative, but can be non-integer).}
  \item{log, log.p}{ logical; if TRUE, probabilities \eqn{p} are given as \eqn{\log(p)}.}
  \item{lower.tail}{ logical; if TRUE (default), probabilities are \eqn{P(X \le x)},
    otherwise, \eqn{P(X > x)}.}
}
\value{
  \code{dchi}, \code{pchi}, and \code{qchi} are respectively the density, distribution 
  function and quantile function of the chi distribution. \code{rchi} generates random 
  deviates drawn from the chi distribution, the length of the result is determined by \code{n}.
}
\details{
  If \code{df} is not specified, they assume the default value of \code{1}.

  The chi distribution with \eqn{\nu} degrees of freedom (also known as shape parameter) has 
  density
  \deqn{
    f(x) =
    \frac{1}{2^{\nu/2-1}\Gamma(\nu/2)} x^{\nu-1}\exp(-x^2/2),
  }
  where \eqn{y > 0}, and \eqn{\nu > 0}. 

  \code{rchi} implements the routine proposed by Monahan (1987) for generating chi-distributed 
  random variables.
}
\author{Felipe Osorio}
\references{
  Forbes, C., Evans, M., Hastings, N., Peacock, B. (2010).
  \emph{Statistical Distributions}, 4th Ed.
  Wiley, New York.

  Johnson, N.L., Kotz, S., Balakrishnan, N. (1994). 
  \emph{Continuous Univariate Distributions, Volume 1}, 2nd Ed. 
  Wiley, New York.

  Monahan, J.F. (1987).
  An algorithm for generating chi random variables.
  \emph{ACM Transactions on Mathematical Software} \bold{13}, 168-172.
}
\seealso{
  \link{Distributions} for other standard distributions.
}
\examples{
x <- rchi(1000, df = 2)
## QQ-plot for chi data against true theoretical distribution:
qqplot(qchi(ppoints(1000), df = 2), x, main = "chi QQ-plot",
  xlab = "Theoretical quantiles", ylab = "Sample quantiles")
abline(c(0,1), lwd = 2, lty = 2, col = "red")
}
\keyword{distribution}
