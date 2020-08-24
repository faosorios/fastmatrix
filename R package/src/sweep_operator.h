/* $ID: sweep_operator.h, last updated 2020-08-16, F.Osorio */

#ifndef FASTMAT_SWEEP_H
#define FASTMAT_SWEEP_H

#include "base.h"

/* sweep operator for symmetric matrices */
extern void sweep_operator(double *, int *, int *, int *, int *, int *);
extern void F77_NAME(sweepop)(double *, int *, int *, int *, int *, int *);

#endif /* FASTMAT_SWEEP_H */
