/* ID: ridge.h, last updated 2020-09-29, F.Osorio */

#ifndef FAST_RIDGE_H
#define FAST_RIDGE_H

/* GCV info required by the minimizer */
typedef struct GCV_info {
  int n, p;
  double edf, pen, GCV, RSS;
  double *u, *d, *y, *rhs, *a, *fitted, *resid;
} GCV_info, *GCVinfo;

#endif /* FAST_RIDGE_H */
