\name{geomean}
\alias{geomean}
\title{Geometric mean}
\description{
  It calculates the geometric mean using a Fused-Multiply-and-Add (FMA) compensated scheme
  for accurate computation of floating-point product.
}
\usage{
geomean(x)
}
\arguments{
  \item{x}{a numeric vector containing the sample observations.}
}
\value{
  The geometric mean of the sample, a non-negative number.
}
\details{
  If \code{x} contains any non-positive values, \code{geomean} returns \code{NA} and
  a warning message is displayed.

  The geometric mean is a measure of central tendency, which is defined as
  \deqn{G = \sqrt[n]{x_1 x_2 \ldots x_n} = \Big(\prod_{i=1}^n x_i\Big)^{1/n}.}

  This procedure calculates the product required in the geometric mean safely using
  a compensated scheme as proposed by Graillat (2009).
}
\seealso{
  \code{\link{mean}}, \code{\link{median}}.
}
\references{
  Graillat, S. (2009).
  Accurate floating-point product and exponentiation.
  \emph{IEEE Transactions on Computers} \bold{58}, 994-1000.

  Oguita, T., Rump, S.M., Oishi, S. (2005).
  Accurate sum and dot product.
  \emph{SIAM Journal on Scientific Computing} \bold{26}, 1955-1988.
}
\examples{
set.seed(149)
x <- rlnorm(1000)
mean(x)    # 1.68169
median(x)  # 0.99663
geomean(x) # 1.01688
}
\keyword{univar}
