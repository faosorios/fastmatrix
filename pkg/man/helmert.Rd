\name{helmert}
\alias{helmert}

\title{Helmert matrix}
\description{
  This function returns the Helmert matrix of order \eqn{n}.
}
\usage{
helmert(n = 1)
}
\arguments{
  \item{n}{ order of the Helmert matrix.}
}
\details{
  A Helmert matrix of order \eqn{n} is a square matrix defined as
  \deqn{\bold{H}_n = \left[{\begin{array}{*{20}{c}}
    1/\sqrt{n} & 1/\sqrt{n} & 1/\sqrt{n} & \dots & 1/\sqrt{n} \\
    1/\sqrt{2} & -1/\sqrt{2} & 0 & \dots & 0 \\
    1/\sqrt{6} & 1/\sqrt{6} & -2/\sqrt{6} & \dots & 0 \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    \frac{1}{\sqrt{n(n-1)}} & \frac{1}{\sqrt{n(n-1)}} & \frac{1}{\sqrt{n(n-1)}} & \dots & -\frac{(n-1)}{\sqrt{n(n-1)}}
  \end{array}}\right].}
  Helmert matrix is orthogonal and is frequently used in the analysis of variance (ANOVA).
}
\value{
  Returns an \eqn{n} by \eqn{n} matrix.
}
\references{
  Lancaster, H.O. (1965).
  The Helmert matrices.
  \emph{The American Mathematical Monthly} \bold{72}, 4-12.

  Gentle, J.E. (2007).
  \emph{Matrix Algebra: Theory, Computations, and Applications in Statistics}.
  Springer, New York.
}
\examples{
n <- 1000
set.seed(149)
x <- rnorm(n)

H <- helmert(n)
object.size(H) # 7.63 Mb of storage
K <- H[2:n,]
z <- c(K \%*\% x)
sum(z^2) # 933.1736

# same that
(n - 1) * var(x)
}
\keyword{array}
\keyword{algebra}
