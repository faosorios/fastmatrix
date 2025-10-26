\name{kronecker.prod}
\alias{kronecker.prod}
\title{Kronecker product on matrices}
\description{
  Computes the kronecker product of two matrices, \eqn{\bold{X}} and \eqn{\bold{Y}}.
}
\usage{
kronecker.prod(x, y = x)
}
\arguments{
  \item{x}{ a numeric matrix or vector.}
  \item{y}{ a numeric matrix or vector.}
}
\details{
  Let \eqn{\bold{X}} be an \eqn{m\times n} and \eqn{\bold{Y}} a \eqn{p\times q} matrix.
  The \eqn{mp\times nq} matrix defined by
  \deqn{\left[{\begin{array}{ccc}
    x_{11}\bold{Y} & \dots & x_{1n}\bold{Y} \\
    \vdots & & \vdots \\
    x_{m1}\bold{Y} & \dots & x_{mn}\bold{Y}
  \end{array}}\right],}
  is called the \emph{Kronecker product} of \eqn{\bold{X}} and \eqn{\bold{Y}}.
}
\value{
  An array with dimensions \code{dim(x) * dim(y)}.
}
\seealso{
  \code{\link{kronecker}} function from \code{base} package is based on \code{\link{outer}}.
  Our \code{C} version is slightly faster.
}
\references{
  Magnus, J.R., Neudecker, H. (2007).
  \emph{Matrix Differential Calculus with Applications in Statistics and Econometrics}, 3rd Edition.
  Wiley, New York.
}
\examples{
# block diagonal matrix:
a <- diag(1:3)
b <- matrix(1:4, ncol = 2)
kronecker.prod(a, b)

# examples with vectors
ones <- rep(1, 4)
y <- 1:3
kronecker.prod(ones, y) # 12-dimensional vector
kronecker.prod(ones, t(y)) # 3 x 3 matrix
}
\keyword{array}
