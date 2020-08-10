# Fast computation of some matrices useful in statistics

Yet another R package for matrices. It contains a small set of functions to fast computation of some matrices and operations useful in statistics.

## Features

Currently we have implemented the following functions:
* Operations envolving the duplication matrix, with minimum requeriments of storage.
* Fast computation of Hadamard product using unrolled loops.
* vec and vech operators to handle rectangular and square matrices.

Our plan in the near future is the implementation of functions to handle:
* Commutation, elimination and symmetrizer matrices.
* Intern products and norms for matrices.
* Array multiplication (see for instance, Appendix A of Wei, 1998).
* Some special matrices arising in numerical analysis.

## Providing Feedback

Please report any bugs/suggestions/improvements to [Felipe Osorio](mailto:felipe.osorios@usm.cl), [Universidad Tecnica Federico Santa Maria](http://www.usm.cl). If you find these routines useful or not then please let me know. Also, acknowledgement of the use of the routines is appreciated.

### To cite the fastmatrix package in publications use:

Osorio, F., Ogueda, A. (2020). fastmatrix: Fast computation of some matrices useful in statistics. 
R package version 0.2. URL: [faosorios.github.io/fastmatrix](https://faosorios.github.io/fastmatrix/)

## About the Authors

Felipe Osorio is an Assistant Professor at [Department of Mathematics](http://www.mat.utfsm.cl/), [Universidad Tecnica Federico Santa Maria](http://www.usm.cl/), Chile.
* Webpage: [fosorios.mat.utfsm.cl](http://fosorios.mat.utfsm.cl/)
* Email: [felipe.osorios [AT] usm.cl](mailto:felipe.osorios@usm.cl)

Alonso Ogueda is a student of the Master of Mathematics offered by the [Department of Mathematics](http://www.mat.utfsm.cl/), [Universidad Tecnica Federico Santa Maria](http://www.usm.cl/), Chile.
* Github: [github.com/aoguedao](https://github.com/aoguedao)


