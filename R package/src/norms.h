/* $ID: norms.h, last updated 2020-08-12, F.Osorio */

#ifndef FASTMAT_NORMS_H
#define FASTMAT_NORMS_H

#include "base.h"

/* vector norms (to be called by R) */
extern void norm_one(double *, int *, int *, double *);
extern void norm_two(double *, int *, int *, double *);
extern void norm_inf(double *, int *, int *, double *);
extern void norm_minkowski(double *, int *, int *, double *, double *);
extern void matrix_norm(double *, int *, int *, int *, int *, double *);

/* to be called by C wrappers */
extern double F77_NAME(minkowski)(int *, double *, int *, double *);
extern double F77_NAME(dnrminf)(int *, double *, int *);

#endif /* FASTMAT_NORMS_H */
