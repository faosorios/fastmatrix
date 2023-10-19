\name{corAR1}
\alias{corAR1}
\title{AR(1) correlation structure}
\usage{
corAR1(rho, p = 2)
}
\arguments{
 \item{rho}{ the value of the lag 1 autocorrelation, which must be between -1 and 1.}
 \item{p}{ dimension of the requested correlation matrix.}
}
\description{
  This function is a constructor for the \code{corAR1} correlation matrix representing
  an autocorrelation structure of order 1.
}
\value{
  Returns a \eqn{p} by \eqn{p} matrix.
}

\examples{
R <- corAR1(rho = 0.8, p = 5)
}
\keyword{array}
