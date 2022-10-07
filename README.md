# Fast computation of some matrices useful in statistics

[![CRAN status](http://www.r-pkg.org/badges/version/fastmatrix)](https://cran.r-project.org/package=fastmatrix)
![CRAN/METACRAN](https://img.shields.io/cran/l/fastmatrix?color=informational)
[![L1pack](https://img.shields.io/badge/Support-L1pack-orange)](https://cran.r-project.org/package=L1pack)
[![SpatialPack](https://img.shields.io/badge/Support-SpatialPack-orange)](https://cran.r-project.org/package=SpatialPack)
![GitHub last commit](https://img.shields.io/github/last-commit/faosorios/fastmatrix)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/fastmatrix)](https://cran.r-project.org/package=fastmatrix)

Yet another R package for matrices. It contains a small set of functions to fast computation of some matrices and operations useful in statistics.

## Features

Latest release (Aug 9, 2022) of [fastmatrix](https://github.com/faosorios/fastmatrix) package have implemented the following functions:
* Array multiplication (see for instance, Appendix A of Wei, 1998).
* Bezier curve based on n+1 control points.
* C version of the Kronecker product which is slightly faster than the built in R base.
* Column-equilibration for rectangular and symmetric matrices.
* Constructors for AR(1) and compound symmetry correlation matrices.
* Covariance matrix estimation using the Mean Square Successive (MSSD) method.
* Estimation of the weighted mean and covariance matrix using an online algorithm (Clarke, 1971).
* Computation of central moments up to fourth order using an online algorithm (Spicer, 1972).
* Fast computation of Hadamard product using unrolled loops.
* Gauss-Seidel, Jacobi and conjugate gradients (CG) iterative methods for solving linear systems.
* Geometric mean using a Fused-Multiply-and-Add (FMA) compensated scheme for accurate computation of floating-point product.
* Inner products and norms for matrices.
* Interface to **C code callable by another C code** from other R packages.
* LDL decomposition for symmetric real matrices.
* Lp norms for vectors.
* LU factorization for square matrices.
* Mahalanobis distances, checking if the covariance is a positive definite matrix.
* Operations envolving the commutation matrix, with minimum requirements of storage.
* Operations envolving the duplication matrix, with minimum requirements of storage.
* Operations envolving the symmetrizer matrix, with minimum requirements of storage.
* Ordinary least-squares (OLS) using several methods: conjugate gradients, Cholesky, QR decomposition, singular value decomposition, and the Sweep operator. This provides an alternative to extend the procedures available in R built-in function 'lm'.
* Power method to compute the dominant eigenvalue and its associated eigenvector.
* Random number generation from the multivariate normal (Gaussian) distribution
* Rank 1 update to Cholesky factorization.
* Ridge estimation for linear regression.
* Routines to compute measures of multivariate skewness and kurtosis proposed by Mardia (1970).
* Routine for the computation of the mediancenter (or geometric median) of multivariate data.
* Routine to compute a Krylov matrix.
* Sherman-Morrison formula.
* Sweep operator for symmetric matrices.
* Test for variance homogeneity of correlated variables (Harris, 1985).
* vec and vech operators to handle rectangular and square matrices.
* Whitening transformation.
* Wilson-Hilferty transformation for chi-squared random variables.

Our plan in the near future is the implementation of functions to handle:
* Some special matrices and operations arising in numerical analysis.

## Reference Manual

* [fastmatrix.pdf](https://cran.r-project.org/web/packages/fastmatrix/fastmatrix.pdf)

## Resources

Latest binaries and sources (development version) for [fastmatrix](https://github.com/faosorios/fastmatrix) are available here (these are local files):

* [fastmatrix_0.4-1245.tar.gz](https://github.com/faosorios/fastmatrix/blob/master/sources/fastmatrix_0.4-1245.tar.gz) - Package sources
* [fastmatrix_0.4-1245.zip](https://github.com/faosorios/fastmatrix/blob/master/binaries/fastmatrix_0.4-1245.zip) - Windows binaries

Version 0.4-1245 of [fastmatrix](https://github.com/faosorios/fastmatrix) can be found at the [CRAN package repository](https://cran.r-project.org/package=fastmatrix):

* [fastmatrix_0.4-1245.tar.gz](https://cran.r-project.org/src/contrib/fastmatrix_0.4-1245.tar.gz) - Package sources
* [fastmatrix_0.4-124.zip](https://cran.r-project.org/bin/windows/contrib/4.2/fastmatrix_0.4-124.zip) - Windows binaries (R-release)
* [fastmatrix_0.4-124.tgz](https://cran.r-project.org/bin/macosx/big-sur-arm64/contrib/4.2/fastmatrix_0.4-124.tgz) - MacOS binaries (R-release, arm64)
* [fastmatrix_0.4-124.tgz](https://cran.r-project.org/bin/macosx/contrib/4.2/fastmatrix_0.4-124.tgz) - MacOS binaries (R-release, x86_64)

## Installation instructions

To install [fastmatrix](https://github.com/faosorios/fastmatrix) **(version 0.4-1245)** from CRAN, start R and enter:
```r
install.packages("fastmatrix")
```

Or install it from its [GitHub repository](https://github.com/faosorios/fastmatrix) **(development version)**. First install the [devtools](https://devtools.r-lib.org/) package.
```r
install.packages("devtools")
```

Then install [fastmatrix](https://github.com/faosorios/fastmatrix) using the `install_github` function in [devtools](https://devtools.r-lib.org/)
```r
library(devtools)
install_github("faosorios/fastmatrix", subdir = "pkg")
```

Alternatively, you can download the source as a tarball (.tar.gz file). Unpack this file (thereby creating a directory named, [fastmatrix](https://github.com/faosorios/fastmatrix)) and install the package source by executing (at the console prompt)
```
R CMD INSTALL fastmatrix
```

Next, you can load the package by using the command: `library(fastmatrix)`

## Providing Feedback

Please report any bugs/suggestions/improvements to [Felipe Osorio](http://fosorios.mat.utfsm.cl/). If you find these routines useful or not then please let me know. Also, acknowledgement of the use of the routines is appreciated.

### To cite the fastmatrix package in publications use:

Osorio, F., Ogueda, A. (2022). fastmatrix: Fast computation of some matrices useful in statistics. 
R package version 0.4-1245. URL: [faosorios.github.io/fastmatrix](https://faosorios.github.io/fastmatrix/)

## About the Authors

Felipe Osorio is an Assistant Professor at [Department of Mathematics](http://www.mat.utfsm.cl/), [Universidad Tecnica Federico Santa Maria](http://www.usm.cl/), Chile.
* Webpage: [fosorios.mat.utfsm.cl](http://fosorios.mat.utfsm.cl/)

Alonso Ogueda is a doctorate student of mathematics from the [Mathematical Sciences Department](https://catalog.gmu.edu/colleges-schools/science/mathematical-sciences/), [George Mason University](https://www2.gmu.edu/), USA.
* Github: [github.com/aoguedao](https://github.com/aoguedao)


