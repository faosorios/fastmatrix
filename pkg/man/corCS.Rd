\name{corCS}
\alias{corCS}
\title{Compound Symmetry Correlation Structure}
\usage{
corCS(rho, p = 2)
}
\arguments{
 \item{rho}{ the value of the correlation between any two correlated observations,
 which must be between -1 and 1.}
 \item{p}{ dimension of the requested correlation matrix.}
}
\description{
  This function is a constructor for the \code{corCS} correlation matrix representing
  a compound symmetry structure corresponding to uniform correlation.
}
\value{
  Returns a \eqn{p} by \eqn{p} matrix.
}

\examples{
R <- corCS(rho = 0.8, p = 5)
}
\keyword{array}
