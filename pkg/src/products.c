/* ID: products.c, last updated 2020-11-03, F.Osorio */

#include "fastmatrix.h"

void
FM_two_product_FMA(double a, double b, double *x, double *y)
{ /* error-free transformation of a product using Fused-Multiply-and-Add
   * Oguita, Rump, Oishi (2005). SIAM Journal on Scientific Computing 26, 1955-1988.
   * doi: 10.1137/030601818 */
  double p;

  *x = p = a * b;
  *y = fma(a, b, -p);
}

void
FM_compensated_product(double *x, int nobs, double *prod)
{ /* product evaluation with a compensated scheme (FMA version)
   * Graillat (2009). IEEE Transactions on Computers 58, 994-1000.
   * doi: 10.1109/TC.2008.215 */
  double accum = 0.0, eps, p, q;

  p = x[0];
  for (int i = 1; i < nobs; i++) {
    FM_two_product_FMA(p, x[i], &q, &eps);
    accum = fma(accum, x[i], eps);
    p = q;
  }
  *prod = p + accum;
}

void
geometric_mean(double *x, int *nobs, double *mean)
{ /* wrapper for 'FM_geometric_mean' */
  int n = *nobs;

  FM_geometric_mean(x, n, mean);
}
