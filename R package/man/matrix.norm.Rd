\name{matrix.norm}
\alias{matrix.norm}
\title{Compute the norm of a rectangular matrix}
\description{
  Computes a matrix norm of \code{x} using LAPACK.  The norm can be the one (\code{"1"})
  norm, the infinity (\code{"inf"}) norm, the Frobenius norm, the maximum modulus
  (\code{"maximum"}) among elements of a matrix, as determined by the value of \code{type}.
}
\usage{
  matrix.norm(x, type = "Frobenius")
}
\arguments{
  \item{x}{ a numeric matrix.}
  \item{type}{ character string, specifying the \emph{type} of matrix norm to be computed.
    A character indicating the type of norm desired.
    \describe{
      \item{\code{"1"}}{ specifies the \bold{one} norm,	(maximum absolute column sum);}
      \item{\code{"Inf"}}{ specifies the \bold{inf}inity norm (maximum absolute row sum);}
      \item{\code{"Frobenius"}}{ specifies the \bold{Frobenius} norm (the	Euclidean norm of
      \code{x} treated as if it were a vector);}
      \item{\code{"maximum"}}{ specifies the \bold{maximum} modulus of all the elements in
      \code{x}.}
    }
}}
\details{
  As function \code{norm} in package \pkg{base}, method of \code{matrix.norm} calls
  the LAPACK function \code{dlange}.

  Note that the 1-, Inf- and maximum norm is faster to calculate than the Frobenius one.
}
\value{
  The matrix norm, a non-negative number.
}
\examples{
# a tiny example
x <- matrix(c(1, 1, 1,
              1, 2, 1,
              1, 3, 1,
              1, 1,-1,
              1, 2,-1,
              1, 3,-1), ncol = 3, byrow = TRUE)
matrix.norm(x, type = "Frobenius")
matrix.norm(x, type = "1")
matrix.norm(x, type = "Inf")

# an example not that small
n <- 1000
x <- .5 * diag(n) + 0.5 * matrix(1, nrow = n, ncol = n)
matrix.norm(x, type = "Frobenius")
matrix.norm(x, type = "1")
matrix.norm(x, type = "Inf")
matrix.norm(x, type = "maximum") # equal to 1
}
\keyword{array}
\keyword{math}
