/* ID: duplication.h, last updated 2020-08-08, F.Osorio */

#ifndef FASTMAT_DUPLICATION_H
#define FASTMAT_DUPLICATION_H

#include "base.h"

/* routines for operations on duplication matrices */
extern void dupl_cols(int *, int *);
extern void duplication_mat(int *, int *, int *, int *);
extern void dupl_left_mult(double *, int *, int *, int *, int *, int *, double *, int *);
extern void dupl_left_trans(double *, int *, int *, int *, int *, int *, int *, double *, int *);
extern void dupl_right_mult(double *, int *, int *, int *, int *, int *, int *, double *, int *);
extern void dupl_right_trans(double *, int *, int *, int *, int *, int *, double *, int *);

#endif /* FASTMAT_DUPLICATION_H */
