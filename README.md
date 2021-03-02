# Fast computation of some matrices useful in statistics

[![CRAN status](http://www.r-pkg.org/badges/version/fastmatrix)](https://cran.r-project.org/package=fastmatrix)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/fastmatrix)](https://cran.r-project.org/package=fastmatrix)

Yet another R package for matrices. It contains a small set of functions to fast computation of some matrices and operations useful in statistics.

## Features

Latest release (February 21, 2021) of **fastmatrix** package have implemented the following functions:
* Array multiplication (see for instance, Appendix A of Wei, 1998).
* C version of the Kronecker product which is slightly faster than the built in R base.
* Column-equilibration for rectangular matrices.
* Covariance matrix estimation using the Mean Square Successive (MSSD) method.
* Estimation of the weighted mean and covariance matrix using an online algorithm (Clarke, 1971).
* Fast computation of Hadamard product using unrolled loops.
* Gauss-Seidel, Jacobi and conjugate gradients (CG) iterative methods for solving linear systems.
* Geometric mean using a Fused-Multiply-and-Add (FMA) compensated scheme for accurate computation of floating-point product.
* Inner products and norms for matrices.
* Interface to **C code callable by another C code** from other R packages.
* LDL decomposition for symmetric real matrices.
* Lp norms for vectors.
* LU factorization for square matrices.
* Mahalanobis distances, checking if the covariance is a positive definite matrix.
* Operations envolving the commutation matrix, with minimum requeriments of storage.
* Operations envolving the duplication matrix, with minimum requeriments of storage.
* Operations envolving the symmetrizer matrix, with minimum requeriments of storage.
* Ordinary least-squares (OLS) using several methods: conjugate gradients, Cholesky, QR decomposition, singular value decomposition, and the Sweep operator. This provides an alternative to extend the procedures available in R built-in function 'lm'.
* Power method to compute the dominant eigenvalue and its associated eigenvector.
* Ridge estimation for linear regression.
* Routines to compute measures of multivariate skewness and kurtosis proposed by Mardia (1970).
* Sherman-Morrison formula.
* Sweep operator for symmetric matrices.
* vec and vech operators to handle rectangular and square matrices.
* Wilson-Hilferty transformation for chi-squared random variables.

Our plan in the near future is the implementation of functions to handle:
* Elimination matrices.
* Some special matrices and operations arising in numerical analysis.

## Reference Manual

<a href="https://github.com/faosorios/fastmatrix/blob/master/man/fastmatrix-manual.pdf"><img src="https://github.com/faosorios/fastmatrix/blob/master/man/fastmatrix.png" height="250"/></a>

## Resources

Latest binaries for **fastmatrix** are available here:

* [fastmatrix_0.3-819.zip](https://github.com/faosorios/fastmatrix/blob/master/binaries/fastmatrix_0.3-819.zip) - Windows binaries
* [fastmatrix_0.3-819.tgz](https://github.com/faosorios/fastmatrix/blob/master/binaries/fastmatrix_0.3-819.tgz) - Mac OS binaries

Version 0.3-81 of **fastmatrix** acan be found at the [CRAN package repository](https://cran.r-project.org/package=fastmatrix)

* [fastmatrix_0.3-81.tar.gz](https://cran.r-project.org/src/contrib/fastmatrix_0.3-81.tar.gz) - Package sources
* [fastmatrix_0.3-81.zip](https://cran.r-project.org/bin/windows/contrib/4.0/fastmatrix_0.3-81.zip) - Windows binaries (R-release)
* [fastmatrix_0.3-81.tgz](https://cran.r-project.org/bin/macosx/contrib/4.0/fastmatrix_0.3-81.tgz) - Mac OS binaries (R-release)

## Installation instructions

To install this package, start R and enter:
```
install.packages("fastmatrix")
```

Alternatively, you can download the source as a tarball or as a zip file. Unpack this file (thereby creating a directory named, fastmatrix) and install the package source by executing (at the console prompt)
```
R CMD INSTALL fastmatrix
```

Next, you can load the package by using the command: `library(fastmatrix)`

## Providing Feedback

Please report any bugs/suggestions/improvements to [Felipe Osorio](http://fosorios.mat.utfsm.cl/). If you find these routines useful or not then please let me know. Also, acknowledgement of the use of the routines is appreciated.

### To cite the fastmatrix package in publications use:

Osorio, F., Ogueda, A. (2021). fastmatrix: Fast computation of some matrices useful in statistics. 
R package version 0.3-819. URL: [faosorios.github.io/fastmatrix](https://faosorios.github.io/fastmatrix/)

## About the Authors

Felipe Osorio is an Assistant Professor at [Department of Mathematics](http://www.mat.utfsm.cl/), [Universidad Tecnica Federico Santa Maria](http://www.usm.cl/), Chile.
* Webpage: [fosorios.mat.utfsm.cl](http://fosorios.mat.utfsm.cl/)

Alonso Ogueda is a student of the Master of Mathematics offered by the [Department of Mathematics](http://www.mat.utfsm.cl/), [Universidad Tecnica Federico Santa Maria](http://www.usm.cl/), Chile.
* Github: [github.com/aoguedao](https://github.com/aoguedao)


