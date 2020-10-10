/* ID: stats.c, last updated 2020-09-03, F.Osorio */

#include "fastmatrix.h"

void
FM_mean_and_var(double *x, int nobs, double *mean, double *var)
{ /* computes the sample mean and variance using an online algorithm */
  int n = 0;
  double accum = 0.0, diff;

  *mean = 0.0;
  for (int i = 0; i < nobs; i++) {
    n++;
    diff = x[i] - *mean;
    *mean += diff / n;
    accum += diff * (x[i] - *mean);
  }
  *var = accum / n;
}

void
FM_online_covariance(double *x, double *y, int nobs, double *xbar, double *ybar,
  double *xvar, double *yvar, double *cov)
{ /* computes the sample covariance using an online algorithm */
  int n = 0;
  double accum = 0.0, acc_x = 0.0, acc_y = 0.0, diff_x, diff_y;

  *xbar = *ybar = 0.0;
  for (int i = 0; i < nobs; i++) {
    n++;
    diff_x = x[i] - *xbar;
    diff_y = y[i] - *ybar;
    *xbar += diff_x / n;
    *ybar += diff_y / n;
    acc_x += diff_x * (x[i] - *xbar);
    acc_y += diff_y * (y[i] - *ybar);
    accum += (n - 1) * (diff_x / n) * (diff_y / n) - accum / n;
  }
  *xvar = acc_x / n;
  *yvar = acc_y / n;
  *cov  = accum;
}

void
FM_center_and_Scatter(double *x, int n, int p, double *weights, double *center, double *Scatter)
{ /* compute center and Scatter estimates using an online algorithm
   * based on AS 41: Applied Statistics 20, 1971, 206-209. doi: 10.2307/2346477 */
  double accum = 0.0, factor = 1.0, wts, *diff;

  /* initialization */
  diff = (double *) Calloc(p, double);
  BLAS1_copy(center, 1, x, n, p); /* copying 1st observation */
  factor /= (double) n;
  BLAS1_scale(factor, weights, 1, n); /* to avoid scaling the Scatter matrix */
  accum += weights[0];

  /* updating stage */
  for (int i = 1; i < n; i++) {
    wts = weights[i];
    accum += wts;
    factor = wts / accum;
    BLAS1_copy(diff, 1, x + i, n, p);
    BLAS1_axpy(-1.0, center, 1, diff, 1, p);
    BLAS1_axpy(factor, diff, 1, center, 1, p);
    factor = wts - factor * wts;
    BLAS2_ger(factor, Scatter, p, p, p, diff, 1, diff, 1);
  }

  Free(diff);
}

void
FM_cov_MSSD(double *x, int n, int p, double *center, double *Scatter)
{ /* compute center and Scatter estimates using the Mean Square Successive Method (MSSD) */
  int accum = 1;
  double *curr, *diff, *prev;

  /* initialization */
  curr = (double *) Calloc(p, double);
  diff = (double *) Calloc(p, double);
  prev = (double *) Calloc(p, double);
  BLAS1_copy(center, 1, x, n, p); /* copying 1st observation */
  BLAS1_copy(prev, 1, x, n, p); /* copying again */

  /* updating stage */
  for (int i = 1; i < n; i++) {
    accum++;
    BLAS1_copy(curr, 1, x + i, n, p); /* current observation */
    Memcpy(diff, curr, p);
    BLAS1_axpy(-1.0, prev, 1, curr, 1, p); /* successive difference */
    BLAS2_ger(0.5 / (n - 1.0), Scatter, p, p, p, curr, 1, curr, 1);
    Memcpy(prev, diff, p);
    BLAS1_axpy(-1.0, center, 1, diff, 1, p);
    BLAS1_axpy(1.0 / accum, diff, 1, center, 1, p);
  }

  Free(curr); Free(diff); Free(prev);
}

double
FM_find_quantile(double *a, int n, int k)
{ /* for an array with n elements, find the element which would be a[k] if
   * the array were sorted from smallest to largest (without the need to do
   * a full sort) */
   double w, x;
   int l = 0;
   int r = n - 1;
   int i, j;

   while (l < r) {
     x = a[k];
     i = l;
     j = r;
     while (j >= i) {
       while (a[i] < x) i++;
       while (x < a[j]) j--;
       if (i <= j) {
         w = a[i];
         a[i++] = a[j];
         a[j--] = w;
       }
     }
     if (j < k) l = i;
     if (k < i) r = j;
  }

  return a[k];
}
