/* ID: commutation.h, last updated 2020-08-25, F.Osorio */

#ifndef FASTMAT_COMMUTATION_H
#define FASTMAT_COMMUTATION_H

#include "base.h"

/* routines for operations on commutation matrices */
extern void F77_NAME(comm_rows)(int *, int *, int *);
extern void F77_NAME(commutation_mat)(int *, int *, int *, int *, int *, int *);
extern void F77_NAME(comm_left_mult)(int *, int *, int *, double *, int *, int *, int *, double *, int *, int *);
extern void F77_NAME(comm_right_mult)(int *, int *, int *, double *, int *, int *, int *, double *, int *, int *);

#endif /* FASTMAT_COMMUTATION_H */
