\name{WH.normal}
\alias{WH.normal}
\title{Wilson-Hilferty transformation for chi-squared variates}
\description{
  Returns the Wilson-Hilferty transformation of random variables with chi-squared distribution.
}
\usage{
WH.normal(x)
}
\arguments{
  \item{x}{ vector or matrix of data with, say, \eqn{p} columns.}
}
\details{
  Let \eqn{T = D^2/p} be a random variable, where \eqn{D^2} denotes the squared Mahalanobis
  distance defined as
  \deqn{D^2 = (\bold{x} - \bold{\mu})^T \bold{\Sigma}^{-1} (\bold{x} - \bold{\mu})}

  Thus the Wilson-Hilferty transformation is given by
  \deqn{z = \frac{T^{1/3} - (1 - \frac{2}{9p})}{(\frac{2}{9p})^{1/2}}}
  and \eqn{z} is approximately distributed as a standard normal distribution. This
  is useful, for instance, in the construction of QQ-plots.
}
\references{
  Wilson, E.B., and Hilferty, M.M. (1931).
  The distribution of chi-square.
  \emph{Proceedings of the National Academy of Sciences of the United States of America} \bold{17}, 684-688.
}
\seealso{\code{\link{Mahalanobis}}}
\examples{
x <- iris[,1:4]
z <- WH.normal(x)
par(pty = "s")
qqnorm(z, main = "Transformed distances QQ-plot")
abline(c(0,1), col = "red", lwd = 2, lty = 2)
}
\keyword{multivariate}
