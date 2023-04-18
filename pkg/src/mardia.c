/* ID: mardia.c, last updated 2020-11-22, F.Osorio */

#include "fastmatrix.h"

void
skewness_and_kurtosis(double *x, int *n, int *p, double *center, double *Scatter, double *stats, int *task)
{ /* wrapper to 'FM_skewness_and_kurtosis' */
  int do_skewness = *task;

  FM_skewness_and_kurtosis(x, *n, *p, center, Scatter, stats, do_skewness);
}
