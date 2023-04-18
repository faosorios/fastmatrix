\name{JarqueBera.test}
\alias{JarqueBera.test}
\title{Jarque-Bera test for univariate normality}
\description{
  Performs an omnibus test for univariate normality.
}
\usage{
JarqueBera.test(x, test = "DH")
}
\arguments{
  \item{x}{ a numeric vector containing the sample observations.}
  \item{test}{ test statistic to be used. One of "DH" (Doornik-Hansen, the default), "JB" (Jarque-Bera), 
  "ALM" (Adjusted Lagrange multiplier).}
}
\value{
  A list of class 'JarqueBera.test' with the following elements:

  \item{statistic}{ value of the statistic, i.e. the value of either Doornik-Hansen, Jarque-Bera, or 
  Adjusted Lagrange multiplier test.}
  \item{parameter}{ the degrees of freedom for the test statistic, which is chi-square distributed.}
  \item{p.value}{ the p-value for the test.}
  \item{skewness}{ the estimated skewness coefficient.}
  \item{kurtosis}{ the estimated kurtosis coefficient.}
  \item{method}{ a character string indicating what type of test was performed.}
}
\references{
  Doornik, J.A., Hansen, H. (2008).
  An omnibus test for univariate and multivariate normality.
  \emph{Oxford Bulletin of Economics and Statistics} \bold{70}, 927-939.

  Jarque, C.M., Bera, A.K. (1980).
  Efficient tests for normality, homoscedasticity and serial independence of regression residuals.
  \emph{Economics Letters} \bold{6}, 255-259.

  Urzua, C.M. (1996).
  On the correct use of omnibus tests for normality.
  \emph{Economics Letters} \bold{53}, 247-251.
}
\examples{
set.seed(149)
x <- rnorm(100)
z <- JarqueBera.test(x, test = "DH")
z

set.seed(173)
x <- runif(100)
z <- JarqueBera.test(x, test = "DH")
z
}
\keyword{htest}
