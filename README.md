# Fast computation of some matrices useful in statistics

[![CRAN status](http://www.r-pkg.org/badges/version/fastmatrix)](https://cran.r-project.org/package=fastmatrix)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/fastmatrix)](https://cran.r-project.org/package=fastmatrix)

Yet another R package for matrices. It contains a small set of functions to fast computation of some matrices and operations useful in statistics.

## Features

Latest release (September 18, 2020) of **fastmatrix** package have implemented the following functions:
* Array multiplication (see for instance, Appendix A of Wei, 1998).
* C version of the Kronecker product which is slightly faster than the built in R base.
* Column-equilibration for rectangular matrices.
* Covariance matrix estimation using the Mean Square Successive (MSSD) method.
* Estimation of the weighted mean and covariance matrix using an online algorithm (Clarke, 1971).
* Fast computation of Hadamard product using unrolled loops.
* Inner products and norms for matrices.
* Interface to **C code callable by another C code** from other R packages.
* Lp norms for vectors.
* Mahalanobis distances, checking if the covariance is a positive definite matrix.
* Operations envolving the commutation matrix, with minimum requeriments of storage.
* Operations envolving the duplication matrix, with minimum requeriments of storage.
* Operations envolving the symmetrizer matrix, with minimum requeriments of storage.
* Ordinary least-squares (OLS) using several methods: Cholesky, QR decomposition, singular value decomposition, and the Sweep operator. This provides an alternative to extend the procedures available in OLS.
* Power method to compute the dominant eigenvalue and its associated eigenvector.
* Procedures for the LU factorization of square matrices.
* Routines to compute measures of multivariate skewness and kurtosis proposed by Mardia (1970).
* Sherman-Morrison formula.
* Sweep operator for symmetric matrices.
* vec and vech operators to handle rectangular and square matrices.

Our plan in the near future is the implementation of functions to handle:
* Elimination matrices.
* Some special matrices and operations arising in numerical analysis.

## Resources

Latest binaries and sources for **fastmatrix** are availables from [CRAN package repository](https://cran.r-project.org/package=fastmatrix)

* [fastmatrix_0.3.tar.gz](https://cran.r-project.org/src/contrib/fastmatrix_0.3.tar.gz) - Package sources
* [fastmatrix_0.3.zip](https://cran.r-project.org/bin/windows/contrib/4.0/fastmatrix_0.3.zip) - Windows binaries (R-release)
* [fastmatrix_0.3.tgz](https://cran.r-project.org/bin/macosx/contrib/4.0/fastmatrix_0.3.tgz) - Mac OS binaries (R-release)
* [fastmatrix.pdf](https://cran.r-project.org/web/packages/fastmatrix/fastmatrix.pdf) - Reference Manual

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

Please report any bugs/suggestions/improvements to [Felipe Osorio](mailto:felipe.osorios@usm.cl), [Universidad Tecnica Federico Santa Maria](http://www.usm.cl). If you find these routines useful or not then please let me know. Also, acknowledgement of the use of the routines is appreciated.

### To cite the fastmatrix package in publications use:

Osorio, F., Ogueda, A. (2020). fastmatrix: Fast computation of some matrices useful in statistics. 
R package version 0.2-3571. URL: [faosorios.github.io/fastmatrix](https://faosorios.github.io/fastmatrix/)

## About the Authors

Felipe Osorio is an Assistant Professor at [Department of Mathematics](http://www.mat.utfsm.cl/), [Universidad Tecnica Federico Santa Maria](http://www.usm.cl/), Chile.
* Webpage: [fosorios.mat.utfsm.cl](http://fosorios.mat.utfsm.cl/)
* Email: [felipe.osorios [AT] usm.cl](mailto:felipe.osorios@usm.cl)

Alonso Ogueda is a student of the Master of Mathematics offered by the [Department of Mathematics](http://www.mat.utfsm.cl/), [Universidad Tecnica Federico Santa Maria](http://www.usm.cl/), Chile.
* Github: [github.com/aoguedao](https://github.com/aoguedao)


