/* ID: mahalanobis.c, last updated 2020-10-17, F.Osorio */

#include "fastmatrix.h"

void
mahal_distances(double *y, int *n, int *p, double *center, double *cov, int *inverted, double *distances)
{ /* wrapper to 'FM_mahal_distances' */
  FM_mahal_distances(y, *n, *p, center, cov, *inverted, distances);
}
