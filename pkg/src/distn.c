/* ID: distn.c, last updated 2025-12-30, F.Osorio */

#include "fastmatrix.h"

/* dpqr-functions for the chi distribution */
static double pdf_chi(double, double, int);
static double cdf_chi(double, double, int, int);
static double quantile_chi(double, double, int, int);
static double rng_chi(double);
/* ..end declarations */

/* ========================================================================== *
 * pdf, cdf, quantile and RNG functions for the chi distribution
 * ========================================================================== */

static double pdf_chi(double x, double df, int log_pdf)
{ /* density of the chi distribution */
  double pdf, val;

  if (x < 0.0) {
    pdf = 0.0;
  } else {
    if (df == 2.0) 
      val = log(x) - 0.5 * SQR(x);
    else {
      val  = (1.0 - 0.5 * df) * M_LN2 - lgammafn(0.5 * df); 
      val += (df - 1.0) * log(x) - 0.5 * SQR(x);
    }
    pdf = exp(val);
  }

  return (log_pdf ? log(pdf) : pdf);
}

static double cdf_chi(double x, double df, int lower, int log_cdf)
{ /* distribution function of the chi distribution */
  double y;

  y = SQR(x);
  return pgamma(y, 0.5 * df, 2.0, lower, log_cdf);
}

static double quantile_chi(double p, double df, int lower, int log_prob)
{ /* quantile function of the chi distribution */
  double y;

  y = qgamma(p, 0.5 * df, 2.0, lower, log_prob);
  return sqrt(y);
}

static double rng_chi(double df)
{ /* random number generation from the chi family of distributions with
   * degrees of freedom parameter df >= 1, using the ratio of uniforms
   * method as described in:
   * Monahan, J.F. (1987). ACM Transactions on Mathematical Software 13, 168-172,
   * Monahan, J.F. (1988). ACM Transactions on Mathematical Software 14, 111. */

  /* Constants : */
  const static double S_EXP_M05    = 0.60653065971263342426; /* exp(-1/2)      */
  const static double S_2_EXP_025  = 2.56805083337548278877; /* 2 * exp(1/4)   */
  const static double S_4_EXP_M135 = 1.03696104258356602834; /* 4 * exp(-1.35) */

  double a, b, c, eta, h, r, u, v, z;

  /* Setup: */
  eta = sqrt(df - 1.0);
  a = S_EXP_M05 * (M_SQRT1_2 + eta) / (0.5 + eta);
  b = -S_EXP_M05 * (1.0 - 0.25 / (SQR(eta) + 1.0));
  c = MAX(-eta, b);

  repeat {
    u = unif_rand();
    v = a + (c - a) * unif_rand();
    z = v / u;

    /* Immediate rejection */
    if (z < -eta)
      continue;

    /* Quick acceptance */
    r = 2.5 - SQR(z);
    if (z < 0.0)
      r += z * SQR(z) / (3.0 * (z + eta));
    if (u < r / S_2_EXP_025)
      return z + eta;

    /* Quick rejection */
    if (SQR(z) > S_4_EXP_M135 / u + 1.4)
      continue;

    /* Regular test */
    h = SQR(eta) * log(1.0 + z / eta) - 0.5 * SQR(z) - z * eta;
    if (2 * log(u) < h)
      return z + eta;
  }
}

void d_chi(int *n, double *y, double *x, double *df, int *ndf, int *give_log)
{ /* pdf of the chi distribution */
  int nobs = *n, nd = *ndf, log_pdf = *give_log;

  for (int i = 0; i < nobs; i++)
    y[i] = pdf_chi(x[i], df[i % nd], log_pdf);
}

void p_chi(int *n, double *y, double *x, double *df, int *ndf, int *lower_tail, int *log_p)
{ /* cdf of the chi distribution */
  int nobs = *n, nd = *ndf, lower = *lower_tail, log_cdf = *log_p;

  for (int i = 0; i < nobs; i++)
    y[i] = cdf_chi(x[i], df[i % nd], lower, log_cdf);
}

void q_chi(int *n, double *y, double *p, double *df, int *ndf, int *lower_tail, int *log_p)
{ /* quantile function of the chi distribution */
  int nobs = *n, nd = *ndf, lower = *lower_tail, log_prob = *log_p;

  for (int i = 0; i < nobs; i++)
    y[i] = quantile_chi(p[i], df[i % nd], lower, log_prob);
}

void r_chi(int *n, double *x, double *df, int *ndf)
{ /* random variates from the chi distribution */
  int nobs = *n, nd = *ndf;

  GetRNGstate();
  for (int i = 0; i < nobs; i++)
    x[i] = rng_chi(df[i % nd]);
  PutRNGstate();
}
