/* $ID: norms.h, last updated 2020-08-11, F.Osorio */

#ifndef FASTMAT_NORMS_H
#define FASTMAT_NORMS_H

#include "base.h"

/* vector norms */
extern void norm_one(double *, int *, int *, double *);
extern void norm_two(double *, int *, int *, double *);
extern void norm_inf(double *, int *, int *, double *);
extern void norm_minkowski(double *, int *, int *, double *, double *);
extern double F77_NAME(minkowski)(int *, double *, int *, double *);

/* matrix norms */
extern void F77_NAME(frobenius_norm)(double *, int *, int *, int *, double *);
extern void F77_NAME(maxcol_norm)(double *, int *, int *, int *, double *);
extern void F77_NAME(maxrow_norm)(double *, int *, int *, int *, double *);

#endif /* FASTMAT_NORMS_H */
