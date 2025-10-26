\name{fastmatrix-package}
\alias{fastmatrix-package}
\alias{fastmatrix.Intro}

\title{\packageTitle{fastmatrix}}

\description{
  \packageDescription{fastmatrix}
}

\details{
  The \pkg{fastmatrix} package provides functions to the efficient construction of \emph{duplication},
  \emph{commutation} and \emph{symmetrizer} matrices with minimal storage requeriments. Common matrix 
  decompositions (e.g. LU, LDL), rank-1 updates (e.g. Cholesky update), iterative solvers for linear 
  systems and other linear algebra utilities, as well as basic matrix operations, such as Hadamard 
  (elementwise) and Kronecker products, the Sherman-Morrison formula and the power method. The package 
  also offers several statistical procedures, namely: the sweep operator, weighted mean and covariance 
  (using online algorithms), ordinary least squares via multiple strategies (Cholesky, QR, SVD, sweep 
  and conjugate gradients), ridge regression (with procedures for selecting the ridge parameter), omnibus
  tests for univariate normality, multivariate skewness and kurtosis measures, Mahalanobis distance (with
  positive-definiteness checks), Wilson-Hilferty transformation of gamma variables (useful, for instance,
  for goodness-of-fit of multivariate normal data), and some random number generators. Finally, the 
  package provides interfaces for code written in C, enabling other R packages (or user-written C code)
  to access the C routines in the \pkg{fastmatrix} package.
}

\author{
  Felipe Osorio \email{faosorios.stat@gmail.com},
  Alonso Ogueda \email{aogueda@gmu.edu}
}

\keyword{package}
