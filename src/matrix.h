/* $ID: matrix.c, last updated 2020-08-08, F.Osorio */

#ifndef FASTMAT_MATRIX_H
#define FASTMAT_MATRIX_H

#include "base.h"

/* operations on vectors */
extern void hadamard_prod(double *, double *, int *, double *);
extern void mat2vech(double *, int *, int *, double *);

#endif /* FASTMAT_MATRIX_H */
