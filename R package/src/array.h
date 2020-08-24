/* $ID: array.h, last updated 2020-08-15, F.Osorio */

#ifndef FASTMAT_ARRAY_H
#define FASTMAT_ARRAY_H

#include "base.h"

/* operations on arrays */
extern void F77_NAME(arraymult)(double *, int *, int *, int *, double *, int *, int *, int *, double *, int *, int *, double *, int *, int *);
extern void F77_NAME(bracketprod)(double *, int *, int *, int *, double *, int *, int *, int *, double *, int *, int *);

#endif /* FASTMAT_ARRAY_H */
