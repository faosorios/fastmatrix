/* $ID: utils.h, last updated 2020-08-16, F.Osorio */

#ifndef FASTMAT_UTILS_H
#define FASTMAT_UTILS_H

#include "base.h"

/* operations on matrices */
extern void F77_NAME(equilibrate_cols)(double *, int *, int *, int *, double *, double *, int *, int *);
extern void F77_NAME(hadamard)(double *, double *, int *, double *);
extern void F77_NAME(inner_frobenius)(double *, int *, double *, int *, int *, int *, double *);

#endif /* FASTMAT_UTILS_H */
