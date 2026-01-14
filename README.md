# fastmatrix: Fast computation of some matrices useful in statistics

[![CRAN status](http://www.r-pkg.org/badges/version/fastmatrix)](https://cran.r-project.org/package=fastmatrix)
![CRAN/METACRAN](https://img.shields.io/cran/l/fastmatrix?color=informational)
![GitHub last commit](https://img.shields.io/github/last-commit/faosorios/fastmatrix)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/fastmatrix)](https://cran.r-project.org/package=fastmatrix)

Yet another R package for matrices. It contains a small set of functions designed to speed up the computation of certain matrix operations that are commonly used in statistics and econometrics. It provides efficient implementations for the computation of several structured matrices, matrix decompositions and statistical procedures, many of which have minimal memory overhead. **Suggested Extension** by [fastverse](https://fastverse.github.io/fastverse/).

Currently, the [fastmatrix](https://github.com/faosorios/fastmatrix) provides support for the following R packages:

[![fChange](https://img.shields.io/badge/Support-fChange-orange)](https://cran.r-project.org/package=fChange)
[![india](https://img.shields.io/badge/Support-india-orange)](https://cran.r-project.org/package=india)
[![L1pack](https://img.shields.io/badge/Support-L1pack-orange)](https://cran.r-project.org/package=L1pack)
[![MVT](https://img.shields.io/badge/Support-MVT-orange)](https://cran.r-project.org/package=MVT)
[![panelSUR](https://img.shields.io/badge/Support-panelSUR-orange)](https://cran.r-project.org/package=panelSUR)
[![RobRegression](https://img.shields.io/badge/Support-RobRegression-orange)](https://cran.r-project.org/package=RobRegression)
[![SAGM](https://img.shields.io/badge/Support-SAGM-orange)](https://cran.r-project.org/package=SAGM)
[![SpatialPack](https://img.shields.io/badge/Support-SpatialPack-orange)](https://cran.r-project.org/package=SpatialPack)

## Resources

Version 0.6-6 of [fastmatrix](https://github.com/faosorios/fastmatrix) can be found at the [CRAN package repository](https://cran.r-project.org/package=fastmatrix):

* [fastmatrix_0.6-6.tar.gz](https://cran.r-project.org/src/contrib/fastmatrix_0.6-6.tar.gz) - Package sources
* [fastmatrix_0.6-4.zip](https://cran.r-project.org/bin/windows/contrib/4.5/fastmatrix_0.6-4.zip) - Windows binaries (R-release)
* [fastmatrix_0.6-6.tgz](https://cran.r-project.org/bin/macosx/big-sur-arm64/contrib/4.5/fastmatrix_0.6-6.tgz) - MacOS binaries (R-release, arm64)
* [fastmatrix_0.6-6.tgz](https://cran.r-project.org/bin/macosx/big-sur-x86_64/contrib/4.5/fastmatrix_0.6-6.tgz) - MacOS binaries (R-release, x86_64)

## Reference Manual

* [PDF manual](https://cran.r-project.org/web/packages/fastmatrix/fastmatrix.pdf) | [HTML manual](https://cran.r-project.org/web/packages/fastmatrix/refman/fastmatrix.html)

## Features

Next, the main functionalities of [fastmatrix](https://github.com/faosorios/fastmatrix) (latest release, version 0.6-6, Jan 14, 2026):
* Structured matrices:
  * Constructors for AR(1) and compound symmetry correlation matrices.
  * Constructors for Frank, Hankel and Helmert matrices.
  * Operations envolving the **commutation matrix**, with minimum requirements of storage.
  * Operations envolving the **duplication matrix**, with minimum requirements of storage.
  * Operations envolving the **symmetrizer matrix**, with minimum requirements of storage.
* Matrix operations and decompositions:
  * Array multiplication (see for instance, Appendix A of [Wei, 1998](https://link.springer.com/book/9789813083295)).
  * C version of the Kronecker product which is slightly faster than the built in R base.
  * Column-equilibration for rectangular and symmetric matrices.
  * Evaluation of a real general matrix polynomial using Horner's scheme.
  * Evaluation of a matrix function where its argument is an upper triangular matrix by applying a Parlett recurrence.
  * Fast computation of Hadamard (elementwise) product using unrolled loops.
  * Gauss-Seidel, Jacobi and conjugate gradients (CG) iterative methods for solving linear systems.
  * Inner products and norms for matrices.
  * Computation of the scaled condition number of a rectangular matrix.
  * LDL decomposition for symmetric real matrices.
  * Computation of the modified Cholesky factorization of a real symmetric but not necessarily positive definite matrix.
  * Lp norms for vectors.
  * LU factorization for square matrices.
  * Matrix square root using the Newton iteration proposed by [Denman and Beavers (1976)](https://doi.org/10.1016/0096-3003(76)90020-5) and the Schur decomposition.
  * Modified Cholesky factorization for symmetric but not necessarily positive definite matrices.
  * Power method to compute the dominant eigenvalue and its associated eigenvector.
  * Rank-1 update to Cholesky factorization.
  * Routine to compute a Krylov matrix.
  * Schur decomposition for a square matrix.
  * Sherman-Morrison formula.
  * Sweep operator for symmetric matrices.
  * vec and vech operators to handle rectangular and square matrices.
* Statistical procedures:
  * Covariance matrix estimation using the Mean Square Successive (MSSD) method.
  * Estimation of the weighted mean and covariance matrix using an online algorithm ([Clarke, 1971](https://doi.org/10.2307/2346477)).
  * Computation of central moments up to fourth order using an online algorithm ([Spicer, 1972](https://doi.org/10.2307/2346477)).
  * Geometric mean using a Fused-Multiply-and-Add (FMA) compensated scheme for accurate computation of floating-point product.
  * Mahalanobis distances, checking if the covariance is a positive definite matrix.
  * Mardia's test for multivariate normality ([Mardia, 1974](https://www.jstor.org/stable/25051892)).
  * Omnibus test for univariate normality ([Jarque-Bera](https://doi.org/10.1016/0165-1765(80)90024-5), [Doornik-Hansen](https://doi.org/10.1111/j.1468-0084.2008.00537.x), Adjusted Lagrange multiplier test ([Urzua, 1996](https://doi.org/10.1016/S0165-1765(96)00923-8)) and a robust version proposed by [Gel and Gastwirt, 2008](https://doi.org/10.1016/j.econlet.2007.05.022)).
  * Ordinary least-squares (OLS) using several methods: conjugate gradients, Cholesky, QR decomposition, singular value decomposition, and the Sweep operator. This provides an alternative to extend the procedures available in R built-in function 'lm'.
  * Ridge estimation for linear regression.
  * Routines to compute measures of multivariate skewness and kurtosis proposed by [Mardia (1970)](https://doi.org/10.2307/2346576).
  * Routine for the computation of the mediancenter (or geometric median) for multivariate data.
  * Test for variance homogeneity of correlated variables ([Harris, 1985](https://doi.org/10.1093/biomet/72.1.103)).
  * Whitening transformation.
  * Wilson-Hilferty transformation for Gamma random variables.
* Random number generators (RNGs):
  * RNG from the multivariate normal (Gaussian) distribution.
  * RNG of uniformly distributed deviates **within** a unitary ball.
  * RNG of uniformly distributed deviats located **on** a spherical surface.
  * RNG of the chi-distribution ([Monahan, 1987](https://doi.org/10.1145/328512.328522))
* Miscellaneous:
  * Bezier curve based on n+1 control points.
  * Floyd-Warshall algorithm to find all shortest paths (if exist) in a directed graph.
* C interfaces:
  * The package provides interfaces for code written in C, **enabling other R packages (or user-written C code) to access the C routines in the fastmatrix package**.

Our plan in the near future is the implementation of functions to handle:
* Some special matrices and operations arising in numerical analysis.

## Installation instructions

To install [fastmatrix](https://github.com/faosorios/fastmatrix) **(version 0.6-6)** from CRAN, start R and enter:
```r
install.packages("fastmatrix")
```

Or install it from its [GitHub repository](https://github.com/faosorios/fastmatrix). First install the [devtools](https://devtools.r-lib.org/) package.
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

Please report any bugs/suggestions/improvements to [Felipe Osorio](https://faosorios.github.io/). If you find these routines useful or not then please let me know. Also, acknowledgement of the use of the routines is appreciated.

### To cite the fastmatrix package in publications use:
``` r
citation("fastmatrix")

To cite fastmatrix in publications use:

  Osorio, F., Ogueda, A. (2026). Fast computation of some matrices
  useful in statistics. R package version 0.6-6. URL:
  https://faosorios.github.io/fastmatrix/

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {Fast computation of some matrices useful in statistics},
    author = {F. Osorio and A. Ogueda},
    year = {2026},
    note = {R package version 0.6-6},
    url = {https://faosorios.github.io/fastmatrix/},
  }
```
## Papers using fastmatrix
- Ogueda, A., Osorio, F. (2025). Influence diagnostics for ridge regression using the Kullback-Leibler divergence. [Statistical Papers](https://doi.org/10.1007/s00362-025-01701-1) 66, 85.
- Guglielmini, S., Claeskens, G. (2025). Asymptotic post-selection inference for regularized graphical models. [Statistics and Computing](https://doi.org/10.1007/s11222-025-10564-3) 35, 36.
- van Buuren, S. (2023). Evaluation and prediction of individual growth trajectories. [Annals of Human Biology](https://doi.org/10.1080/03014460.2023.2190619) 50, 247–257.

## About the Authors

Felipe Osorio is an applied statistician and creator of several R packages
* Webpage: [fosorios.github.io](https://faosorios.github.io/)

Alonso Ogueda is a doctorate student of mathematics from the [Mathematical Sciences Department](https://catalog.gmu.edu/colleges-schools/science/mathematical-sciences/), [George Mason University](https://www2.gmu.edu/), USA.
* Github: [github.com/aoguedao](https://github.com/aoguedao)
