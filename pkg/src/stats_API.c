/* ID: stats_API.c, last updated 2023-07-23, F.Osorio */

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
FM_moments(double *x, int nobs, double *mean, double *s2, double *s3, double *s4)
{ /* compute sample mean and sums of powers of deviations, this is a slight
   * modification of AS 52: Applied Statistics 21, 1972, 226-227.
   * doi: 10.2307/2346507 */
  double accum1, accum2, accum3, accum4, m1, m2, n = 1.0, diff, term;

  /* initialization */
  accum1 = x[0];
  accum2 = accum3 = accum4 = 0.0;

  /* updating stage */
  for (int i = 1; i < nobs; i++) {
    n++;
    diff = x[i] - accum1;
    term = diff / n;
    m1 = n - 1.;
    m2 = n - 2.;
    accum4 -= term * (4. * accum3 - term * (6. * accum2 + m1 * (1. + CUBE(m1)) * SQR(term)));
    accum3 -= term * (3. * accum2 - n * m1 * m2 * SQR(term));
    accum2 += n * m1 * SQR(term);
    accum1 += term;
  }

  /* saving results */
  *mean = accum1;
  *s2   = accum2 / n;
  *s3   = accum3 / n;
  *s4   = accum4 / n;
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
FM_geometric_mean(double *x, int nobs, double *mean)
{ /* computes the geometric mean using a compensated product scheme */
  FM_compensated_product(x, nobs, mean);
}

void
FM_online_center(double *x, int n, int p, double *weights, double *center)
{ /* compute center estimate using an online algorithm
   * based on AS 41: Applied Statistics 20, 1971, 206-209.
   * doi: 10.2307/2346477 */
  double accum = 0.0, factor = 1.0, wts, *diff, *mean;

  /* initialization */
  diff = (double *) Calloc(p, double);
  mean = (double *) Calloc(p, double);
  BLAS1_copy(mean, 1, x, n, p);
  accum += weights[0];

  /* updating stage */
  for (int i = 1; i < n; i++) {
    wts = weights[i];
    accum += wts;
    factor = wts / accum;
    BLAS1_copy(diff, 1, x + i, n, p);
    BLAS1_axpy(-1.0, mean, 1, diff, 1, p);
    BLAS1_axpy(factor, diff, 1, mean, 1, p);
  }

  /* saving results */
  BLAS1_copy(center, 1, mean, 1, p);

  Free(diff); Free(mean);
}

void
FM_center_and_Scatter(double *x, int n, int p, double *weights, double *center, double *Scatter)
{ /* compute center and Scatter estimates using an online algorithm
   * based on AS 41: Applied Statistics 20, 1971, 206-209.
   * doi: 10.2307/2346477 */
  double accum = 0.0, factor = 1.0, wts, *diff, *mean, *cov;

  /* initialization */
  diff = (double *) Calloc(p, double);
  mean = (double *) Calloc(p, double);
  cov  = (double *) Calloc(p * p, double);
  BLAS1_copy(mean, 1, x, n, p); /* copying 1st observation */
  accum += weights[0];

  /* updating stage */
  for (int i = 1; i < n; i++) {
    wts = weights[i];
    accum += wts;
    factor = wts / accum;
    BLAS1_copy(diff, 1, x + i, n, p);
    BLAS1_axpy(-1.0, mean, 1, diff, 1, p);
    BLAS1_axpy(factor, diff, 1, mean, 1, p);
    factor = wts - factor * wts;
    BLAS2_ger(factor, cov, p, p, p, diff, 1, diff, 1);
  }

  /* saving results */
  BLAS1_copy(center, 1, mean, 1, p);
  factor = 1.0 / (double) n;
  FM_scale_mat(Scatter, p, factor, cov, p, p, p);

  Free(diff); Free(mean); Free(cov);
}

void
FM_skewness_and_kurtosis(double *x, int n, int p, double *center, double *Scatter, double *stats, int do_skewness)
{ /* computes Mardia's multivariate skewness and kurtosis */
  char *side = "R", *uplo = "L", *trans = "T", *diag = "N";
  int info = 0, job = 0;
  double dist, skew = 0.0, kurt = 0.0;

  /* computes the triangular factor of 'Scatter' matrix */
  FM_chol_decomp(Scatter, p, p, job, &info);
  if (info)
    error("Covariance matrix is possibly not positive-definite");

  /* standardizing the rows of the data matrix */
  FM_centering(x, n, p, center);
  BLAS3_trsm(1.0, Scatter, p, n, p, side, uplo, trans, diag, x, n);

  /* computation of kurtosis coefficient */
  for (int i = 0; i < n; i++) {
    dist  = FM_norm_sqr(x + i, n, p);
    skew += R_pow_di(dist, 3);
    kurt += SQR(dist);
  }

  if (!do_skewness) { /* skewness coefficient is not required */
    stats[0] = 0.0;
    stats[1] = kurt / n;
    return;
  }

  /* computation of skewness coefficient */
  for (int i = 0; i < n; i++){
    for (int j = i + 1; j < n; j++) {
      dist  = BLAS1_dot_product(x + i, n, x + j, n, p);
      skew += 2.0 * R_pow_di(dist, 3);
    }
  }

  /* copying coefficients */
  stats[0] = skew / SQR(n);
  stats[1] = kurt / n;
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

void
FM_cov4th(double *x, int n, int p, double *center, double *cov)
{ /* distribution-robust (correction factor of the) covariance matrix
   * (see Section 5 of Harris (1985). Biometrika 72, 103-107) */
  double accum, dr, ds;

  for (int r = 0; r < p; r++) {
    for (int s = r; s < p; s++) {
      accum = 0.0;
      for (int i = 0; i < n; i++) {
        dr = x[i + r * n] - center[r];
        ds = x[i + s * n] - center[s];
        accum += SQR(dr) * SQR(ds);
      }
      *(cov + r + s * p) = accum / n;
      *(cov + s + r * p) = *(cov + r + s * p);
    }
  }
}
