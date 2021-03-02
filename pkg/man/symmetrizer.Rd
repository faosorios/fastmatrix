\name{symmetrizer}
\alias{symmetrizer}

\title{Symmetrizer matrix}
\description{
  This function returns the symmetrizer matrix of order \eqn{n} which transforms,
  for every \eqn{n\times n} matrix \eqn{\bold{A}}, \code{vec}\eqn{(\bold{A})} into
  \code{vec}\eqn{((\bold{A} + \bold{A}^T)/2)}.
}
\usage{
symmetrizer(n = 1, matrix = FALSE)
}
\arguments{
  \item{n}{ order of the symmetrizer matrix.}
  \item{matrix}{ a logical indicating whether the symmetrizer matrix will be returned.}
}
\details{
  This function is a wrapper function for the function \code{symm.info}. This function
  provides the information required to create the symmetrizer matrix. If option \code{matrix = FALSE}
  the symmetrizer matrix is stored in three vectors containing the coordinate list of
  indexes for rows, columns and the values.

  \bold{Warning:} \code{matrix = TRUE} is \bold{not} recommended, unless the order \code{n}
  be small. This matrix can require a huge amount of storage.
}
\value{
  Returns an \eqn{n^2} by \eqn{n^2} matrix (if requested).
}
\references{
  Abadir, K.M., Magnus, J.R. (2005).
  \emph{Matrix Algebra}.
  Cambridge University Press.

  Magnus, J.R., Neudecker, H. (2007).
  \emph{Matrix Differential Calculus with Applications in Statistics and Econometrics}, 3rd Edition.
  Wiley, New York.
}
\seealso{
  \code{\link{symm.info}}
}
\examples{
z <- symmetrizer(n = 100)
object.size(z) # 319 Kb of storage

N100 <- symmetrizer(n = 100, matrix = TRUE) # time: < 2 secs
object.size(N100) # 800 Mb of storage, do not request this matrix!

# a small example
N3 <- symmetrizer(n = 3, matrix = TRUE)
a <- matrix(rep(c(2,4,6), each = 3), ncol = 3)
a
b <- 0.5 * (a + t(a))
b
v <- N3 \%*\% vec(a)
all(vec(b) == as.vector(v)) # vectors are equal!
}
\keyword{array}
\keyword{algebra}
