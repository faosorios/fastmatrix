/* $ID: matrix.h, last updated 2020-08-08, F.Osorio */

#ifndef FASTMAT_MATRIX_H
#define FASTMAT_MATRIX_H

#include "base.h"

/* operations on matrices */
extern void equilibrate(double *, int *, int *, int *, double *, double *, int *);
extern void mat2vech(double *, int *, int *, double *);

#endif /* FASTMAT_MATRIX_H */
