/* ID: base.h, last updated 2020-08-01, F.Osorio */

#ifndef FASTMAT_BASE_H
#define FASTMAT_BASE_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Print.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

/* some definitions */
#define DNULLP   (double *) 0
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define SQR(x)   R_pow_di(x, 2)
#define repeat   for(;;)

#endif /* FASTMAT_BASE_H */
