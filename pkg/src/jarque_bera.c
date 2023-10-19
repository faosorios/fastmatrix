/* ID: jarque_bera.c, last updated 2023-07-23, F.Osorio */

#include "fastmatrix.h"

/* static functions */
static double skewness_to_z1(double, double);
static double kurtosis_to_z2(double, double, double);
static double sqr_b1_ALM(double, double);
static double sqr_b2_ALM(double, double);
/* ..end declarations */

void
doornik_hansen(double *x, int *nobs, double *skew, double *kurt, double *stat)
{ /* computes an omnibus test for univariate normality.
   * Doornik and Hansen. Oxford B. Econ. Stat. 70, 927-939, 2008.
   * doi: 10.1111/j.1468-0084.2008.00537.x */
  int n = *nobs;
  double mean, m2, m3, m4, b1, b2, z1, z2;

  /* computes central moments */
  FM_moments(x, n, &mean, &m2, &m3, &m4);

  /* skewness and kurtosis */
  b1 = m3 / R_pow(m2, 1.5);
  b2 = m4 / SQR(m2);

  /* applying transformations to achieve normality */
  z1 = skewness_to_z1(b1, (double) n);
  z2 = kurtosis_to_z2(b1, b2, (double) n);

  /* saving statistics */
  *skew = b1;
  *kurt = b2;

  /* Doornik-Hansen test statistic */
  *stat = SQR(z1) + SQR(z2);
}

double 
skewness_to_z1(double skew, double n)
{ /* transform the skewness into a standard normal using the 
   * Johnson SU approximation, as proposed by D'Agostino
   * Biometrika 57, 679-681, 1970. */
  double a, b, w2, d, y, z;

  b  = 3. * (SQR(n) + 27. * n - 70.) * (n + 1.) * (n + 3.) /
       ((n - 2.) * (n + 5.) * (n + 7.) * (n + 9.));
  w2 = -1. + sqrt(2. * (b - 1.));
  d  = 1. / sqrt(log(sqrt(w2)));
  y  = skew * sqrt((n + 1.) * (n + 3.) / (6. * (n - 2.)));
  a  = sqrt(2. / (w2 - 1.));
  z  = d * log((y / a) + sqrt(SQR(y / a) + 1.));

  return z;
}

double 
kurtosis_to_z2(double skew, double kurt, double n)
{ /* transform the kurtosis into a standard normal using the 
   * Wilson-Hilferty cubed root transformation */
  double a, d, c, k, alpha, p1, p2, p3, p4, q, z;
  double n2 = SQR(n), n3 = CUBE(n);

  p1 = n2 + 15. * n - 4.;
  p2 = n2 + 27. * n - 70.;
  p3 = n2 + 2. * n - 5.;
  p4 = n3 + 37. * n2 + 11. * n - 313.;

  d = (n - 3.) * (n + 1.) * p1;
  a = (n - 2.) * (n + 5.) * (n + 7.) * p2 / (6. * d);
  c = (n - 7.) * (n + 5.) * (n + 7.) * p3 / (6. * d);
  k = (n + 5.) * (n + 7.) * p4 / (12. * d);

  alpha = a + SQR(skew) * c;
  q = 2. * (kurt - 1. - SQR(skew)) * k;

  z  = R_pow(0.5 * q / alpha, 1. / 3.) - 1. + 1. / (9. * alpha);
  z *= sqrt(9. * alpha);

  return z;
}

void
jarque_bera(double *x, int *nobs, double *skew, double *kurt, double *stat)
{ /* computes Jarque-Bera test for normality based on the standardized 
   * third and fourth moments.
   * Jarque and Bera. Economics Letters 6, 255-259, 1980. 
   * doi: 10.1016/0165-1765(80)90024-5 */
  int n = *nobs;
  double mean, m2, m3, m4, b1, b2;

  /* computes central moments */
  FM_moments(x, n, &mean, &m2, &m3, &m4);

  /* skewness and kurtosis */
  b1 = m3 / R_pow(m2, 1.5);
  b2 = m4 / SQR(m2);

  /* saving statistics */
  *skew = b1;
  *kurt = b2;

  /* Jarque-Bera test statistic */
  *stat = ((double) n / 6.) * (SQR(b1) + SQR(b2 - 3.) / 4.);
}

void
robust_JB(double *x, double *z, int *nobs, double *skew, double *kurt, double *stat)
{ /* computes a robust Jarque-Bera test for normality based on the average 
   * absolute deviation from the sample median (MAAD).
   * Gel and Gastwirth. Economics Letters 99, 30-32, 2008. 
   * doi: 10.1016/j.econlet.2007.05.022 */
  int n = *nobs;
  double mad, mean, m2, m3, m4, b1, b2, c1 = 6., c2 = 64.;

  /* computes central moments and MAAD */
  FM_moments(x, n, &mean, &m2, &m3, &m4);
  mad = BLAS1_sum_abs(z, 1, n) / n;
  mad /= M_SQRT_2dPI; /* MAAD correction */

  /* skewness and kurtosis */
  b1 = m3 / CUBE(mad);
  b2 = m4 / FOURTH(mad);

  /* saving statistics */
  *skew = b1;
  *kurt = b2;

  /* Jarque-Bera test statistic */
  *stat = ((double) n) * (SQR(b1) / c1 + SQR(b2 - 3.) / c2);
}

void
urzua_ALM(double *x, int *nobs, double *skew, double *kurt, double *stat)
{ /* computes the adjusted Lagrange multiplier test for normality.
   * Urzua. Economics Letters 53, 247-251, 1996. 
   * doi: 10.1016/S0165-1765(96)00923-8 */
  int n = *nobs;
  double mean, m2, m3, m4, b1, b2, sqr1, sqr2;

  /* computes central moments */
  FM_moments(x, n, &mean, &m2, &m3, &m4);

  /* skewness and kurtosis */
  b1 = m3 / R_pow(m2, 1.5);
  b2 = m4 / SQR(m2);

  /* applying transformations to achieve chi-squared distributions */
  sqr1 = sqr_b1_ALM(b1, (double) n);
  sqr2 = sqr_b2_ALM(b2, (double) n);

  /* saving statistics */
  *skew = b1;
  *kurt = b2;

  /* adjusted Lagrange multiplier test statistic */
  *stat = sqr1 + sqr2;
}

double 
sqr_b1_ALM(double skew, double n)
{ /* exact variance of the skewness under normality */
  double sqr, var;

  var = 6. * (n - 2.) / ((n + 1.) * (n + 3.));
  sqr = SQR(skew) / var; 

  return sqr;
}

double 
sqr_b2_ALM(double kurt, double n)
{ /* exact mean and variance of kurtosis under normality */
  double sqr, mean, var;

  mean = 3. * (n - 1.) / (n + 1.);
  var  = 24. * n * (n - 2.) * (n - 3.) / (SQR(n + 1.) * (n + 3.) * (n + 5.));
  sqr  = SQR(kurt - mean) / var; 

  return sqr;
}
