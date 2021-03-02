\name{commutation}
\alias{commutation}
\title{Commutation matrix}
\description{
  This function returns the commutation matrix of order \eqn{mn} which transforms,
  for an \eqn{m\times n}  matrix \eqn{\bold{A}}, \code{vec}\eqn{(\bold{A})} to
  \code{vec}\eqn{(\bold{A}^T)}.
}
\usage{
commutation(m = 1, n = m, matrix = FALSE, condensed = FALSE)
}
\arguments{
  \item{m}{ a positive integer row dimension.}
  \item{n}{ a positive integer column dimension.}
  \item{matrix}{ a logical indicating whether the commutation matrix will be returned.}
  \item{condensed}{ logical. Information should be returned in compact form?}
}
\details{
  This function is a wrapper function for the function \code{comm.info}. This function
  provides the minimum information required to create the commutation matrix. If option
  \code{matrix = FALSE} the commutation matrix is stored in two vectors containing the
  coordinate list of indexes for rows and columns. Option \code{condensed = TRUE} only
  returns vector of indexes for the rows of commutation matrix.

  \bold{Warning:} \code{matrix = TRUE} is \bold{not} recommended, unless the order \code{m}
  \bold{and} \code{n} be small. This matrix can require a huge amount of storage.
}
\value{
  Returns an \eqn{mn} by \eqn{mn} matrix (if requested).
}
\references{
  Magnus, J.R., Neudecker, H. (1979).
  The commutation matrix: some properties and applications.
  \emph{The Annals of Statistics} \bold{7}, 381-394.

  Magnus, J.R., Neudecker, H. (2007).
  \emph{Matrix Differential Calculus with Applications in Statistics and Econometrics}, 3rd Edition.
  Wiley, New York.
}
\seealso{
  \code{\link{comm.info}}
}
\examples{
z <- commutation(m = 100, condensed = TRUE)
object.size(z) # 40.6 Kb of storage

z <- commutation(m = 100, condensed = FALSE)
object.size(z) # 80.7 Kb of storage

K100 <- commutation(m = 100, matrix = TRUE) # time: < 2 secs
object.size(K100) # 400 Mb of storage, do not request this matrix!

# a small example
K32 <- commutation(m = 3, n = 2, matrix = TRUE)
a <- matrix(1:6, ncol = 2)
v <- K32 \%*\% vec(a)
all(vec(t(a)) == as.vector(v)) # vectors are equal!
}
\keyword{array}
\keyword{algebra}
