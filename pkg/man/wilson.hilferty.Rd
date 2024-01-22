\name{wilson.hilferty}
\alias{wilson.hilferty}
\title{Wilson-Hilferty transformation}
\description{
  Returns the Wilson-Hilferty transformation of random variables with Gamma distribution.
}
\usage{
  wilson.hilferty(x, shape, rate = 1)
}
\arguments{
  \item{x}{ a numeric vector containing Gamma distributed deviates.}
  \item{shape, rate}{shape and rate parameters. Must be positive.}
}
\details{
  Let \eqn{X} be a random variable following a Gamma distribution with parameters \eqn{a} = \code{shape} 
  and \eqn{b} = \code{rate} with density
  \deqn{
    f(x) = \frac{b^a}{\Gamma(a)} x^{a-1}\exp(-bx),}{
    f(x) = b^a/Gamma(a) x^(a-1) exp(-bx),
  }
  where \eqn{x \ge 0}, \eqn{a > 0}{a > 0}, \eqn{b > 0}{b > 0} and consider the random variable 
  \eqn{T = X/(a/b)}. Thus, the Wilson-Hilferty transformation
  \deqn{z = \frac{T^{1/3} - (1 - \frac{1}{9a})}{(\frac{1}{9a})^{1/2}}}
  is approximately distributed as a standard normal distribution. This is useful, for instance, 
  in the construction of QQ-plots.
}
\references{
  Terrell, G.R. (2003).
  The Wilson-Hilferty transformation is locally saddlepoint.
  \emph{Biometrika} \bold{90}, 445-453.

  Wilson, E.B., and Hilferty, M.M. (1931).
  The distribution of chi-square.
  \emph{Proceedings of the National Academy of Sciences of the United States of America} \bold{17}, 684-688.
}
\seealso{\code{\link{WH.normal}}}
\examples{
x <- rgamma(n = 300, shape = 2, rate = 1)
z <- wilson.hilferty(x, shape = 2, rate = 1)
par(pty = "s")
qqnorm(z, main = "Transformed Gamma QQ-plot")
abline(c(0,1), col = "red", lwd = 2, lty = 2)
}
\keyword{univar}
