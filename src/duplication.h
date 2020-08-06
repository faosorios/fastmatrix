/* ID: duplication.h, last updated 2020-08-03, F.Osorio */

#ifndef FASTMAT_DUPLICATION_H
#define FASTMAT_DUPLICATION_H

#include "base.h"

extern void duplication_mat(int *, int *, int *, int *);
extern void dupl_left_mult(double *, int *, int *, int *, int *, int *, double *, int *);
extern void dupl_left_trans(double *, int *, int *, int *, int *, int *, int *, double *, int *);
extern void dupl_right_trans(double *, int *, int *, int *, int *, int *, double *, int *);

#endif /* FASTMAT_DUPLICATION_H */
