\name{duplication}
\alias{duplication}

\title{Duplication matrix}
\description{
  This function returns a duplication matrix of order \code{n} which transforms,
  for a symmetric matrix \code{x}, \code{vech(x)} into \code{vec(x)}.
}
\usage{
duplication(n = 1)
}
\arguments{
  \item{n}{ order of the duplication matrix.}
}
\value{
  Returns an \eqn{n^2} by \eqn{n(n + 1)/2} matrix.
}
\references{
  Magnus, J.R., and Neudecker, H. (1980).
  The elimination matrix, some lemmas and applications.
  \emph{SIAM Journal on Algebraic Discrete Methods} \bold{1}, 422-449.

  Magnus, J.R., and Neudecker, H. (1999).
  \emph{Matrix Differential Calculus with Applications in Statistics and Econometrics}, 2nd Edition.
  Wiley, New York.
}
\examples{
D3 <- duplication(3, matrix = TRUE)
a <- matrix(c( 1, 2, 3,
               2, 3, 4,
               3, 4, 5), nrow = 3)
vh <- a[upper.tri(a, diag = TRUE)]
v <- D3 \%*\% vh
as.vector(vh)
}
\keyword{array}
\keyword{algebra}
