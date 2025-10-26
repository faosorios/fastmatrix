\name{mardia.test}
\alias{mardia.test}
\title{Mardia test for multivariate normality}
\description{
  Performs Mardia's tests to assess multivariate normality based on the multivariate skewness and kurtosis coefficients.
}
\usage{
mardia.test(x)
}
\arguments{
  \item{x}{ matrix of data with, say, \eqn{p} columns.}
}
\value{
  A list of class 'Mardia.test' with the following elements:

  \item{skewness}{ a list containig the \code{statistic}, the degrees of freedom (\code{parameter}), 
    the \code{p.value} for the test and the estimated \code{skewness} coefficient.}
  \item{kurtosis}{ a list containig the \code{statistic}, the \code{p.value} for the test and the 
    estimated \code{kurtosis} coefficient.}
}
\references{
  Mardia, K.V. (1970).
  Measures of multivariate skewness and kurtosis with applications.
  \emph{Biometrika} \bold{57}, 519-530.

  Mardia, K.V. (1974).
  Applications of some measures of multivariate skewness and kurtosis in testing normality and robustness studies.
  \emph{Sankhya} \bold{36}, 115-128.
}
\examples{
setosa <- iris[1:50,1:4]
z <- mardia.test(setosa)
z

set.seed(149)
Sigma <- matrix(c(10,3,3,2), ncol = 2)
x <- rmnorm(n = 300, Sigma = Sigma)
z <- mardia.test(x)
z
}
\keyword{htest}
\keyword{multivariate}
