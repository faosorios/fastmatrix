\name{frank}
\alias{frank}

\title{Frank matrix}
\description{
  This function returns the Frank matrix of order \eqn{n}.
}
\usage{
frank(n = 1)
}
\arguments{
  \item{n}{ order of the Frank matrix.}
}
\details{
  A Frank matrix of order \eqn{n} is a square matrix \eqn{\bold{F}_n = (f_{ij})} defined as
  \deqn{f_{ij} = \left\{ {\begin{array}{ll}
    n - j + 1, & i \le j, \\
    n - j, & i = j + 1, \\
    0, & i \ge j + 2
  \end{array}} \right.}
}
\value{
  Returns an \eqn{n} by \eqn{n} matrix.
}
\references{
  Frank, W.L. (1958).
  Computing eigenvalues of complex matrices by determinant evaluation and by methods of Danilewski and Wielandt.
  \emph{Journal of the Society for Industrial and Applied Mathematics} \bold{6}, 378-392.
}
\examples{
F5 <- frank(n = 5)
det(F5) # equals 1
}
\keyword{array}
\keyword{algebra}
