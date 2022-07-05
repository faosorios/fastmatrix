\name{harris.test}
\alias{harris.test}
\title{Test for variance homogeneity of correlated variables}
\description{
  Performs large-sample methods for testing equality of \eqn{p \ge 2} correlated variables.
}
\usage{
harris.test(x, test = "Wald")
}
\arguments{
  \item{x}{a matrix or data frame. As usual, rows are observations and columns are variables.}
  \item{test}{ test statistic to be used. One of "Wald" (default), "log", "robust" or "log-robust".}
}
\value{
  A list of class 'harris.test' with the following elements:

  \item{statistic}{ value of the statistic, i.e. the value of either Wald test, using the log-transformation,
  or distribution-robust versions of the test (robust and log-robust).}
  \item{parameter}{ the degrees of freedom for the test statistic, which is chi-square distributed.}
  \item{p.value}{ the p-value for the test.}
  \item{estimate}{ the estimated covariance matrix.}
  \item{method}{ a character string indicating what type of test was performed.}
}
\references{
  Harris, P. (1985).
  Testing the variance homogeneity of correlated variables.
  \emph{Biometrika} \bold{72}, 103-107.
}
\examples{
x <- iris[,1:4]
z <- harris.test(x, test = "robust")
z
}
\keyword{htest}
