/* $ID: utils.h, last updated 2020-08-16, F.Osorio */

#ifndef FASTMAT_UTILS_H
#define FASTMAT_UTILS_H

#include "base.h"

/* operations on vectors */
extern void F77_NAME(hadamard)(double *, double *, int *, double *);
extern void F77_NAME(inner_frobenius)(double *, int *, double *, int *, int *, int *, double *);
extern void F77_NAME(sweepop)(double *, int *, int *, int *, double *, int *, int *);
extern void F77_NAME(vech)(double *, int *, int *, double *);

#endif /* FASTMAT_UTILS_H */
