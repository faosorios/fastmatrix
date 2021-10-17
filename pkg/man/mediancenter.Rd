\name{mediancenter}
\alias{mediancenter}
\title{Mediancenter}
\description{
  It calculates the mediancenter (or geometric median) of multivariate data.
}
\usage{
mediancenter(x)
}
\arguments{
  \item{x}{a matrix or data frame. As usual, rows are observations and columns are
    variables.}
}
\value{
  A list containing the following named components:
  \item{median}{an estimate for the mediancenter of the data.}
  \item{iter}{the number of iterations performed, it is negative if a degenerate
              solution is found.}
}
\details{
  The mediancenter for a sample of multivariate observations is computed using a
  steepest descend method combined with bisection. The mediancenter invariant to
  rotations of axes and is useful as a multivariate generalization of the median
  of univariate sample.
}
\seealso{
  \code{\link{cov.wt}}, \code{\link{median}}.
}
\references{
  Gower, J.C. (1974).
  Algorithm AS 78: The mediancentre.
  \emph{Applied Statistics} \bold{23}, 466-470.
}
\examples{
x <- cbind(1:10, c(1:3, 8:5, 8:10))
z <- mediancenter(x)$median # degenerate solution
xbar <- colMeans(x)
plot(x, xlab = "", ylab = "")
points(x = xbar[1], y = xbar[2], pch = 16, col = "red")
points(x = z[1], y = z[2], pch = 3, col = "blue", lwd = 2)
}
\keyword{multivariate}
