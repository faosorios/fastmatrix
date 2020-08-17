/* $ID: matrix.h, last updated 2020-08-08, F.Osorio */

#ifndef FASTMAT_MATRIX_H
#define FASTMAT_MATRIX_H

#include "base.h"

/* operations on vectors */
extern void hadamard_prod(double *, double *, int *, double *);
extern void mat2vech(double *, int *, int *, double *);
extern void power_method(double *, int *, int *, int *, double *, double *, int *, double *, int *);

#endif /* FASTMAT_MATRIX_H */
