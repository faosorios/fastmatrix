/* $ID: brent.c, last updated 2020-11-26, F.Osorio */

#include "fastmatrix.h"

double FM_brent(double ax, double bx, double (*f)(double, void *), void *info, double tolerance)
{ /* Brent's method for unidimensional minimization.
   * This is a recoded version of fmin.f extracted from netlib.org which is
   * a slightly modified version of the 'localmin' procedure described in
   * Brent, R.P. (1973). Algorithms for Minimization Without Derivatives.
   * Dover, New York, pp. 79-80 and 188-190. */
  double a = ax, b = bx;
  double d, e, m, p, q, r, t2, u, v, w, x, fu, fv, fw, fx;
  double eps = DOUBLE_EPS; /* machine epsilon */
  double tol1 = eps + 1.0, tol3 = tolerance / 3.0;

  /* initialization */
  eps  = sqrt(eps);
  d = e = 0.0;
  v = w = x = a + GOLDEN * (b - a);
  fx = (*f)(x, info);
  fv = fw = fx;

  /* main loop */
  repeat {
    m = 0.5 * (a + b);
    tol1 = eps * fabs(x) + tol3;
    t2 = 2.0 * tol1;

    /* check stopping criterion */
    if (fabs(x - m) <= t2 - 0.5 * (b - a))
      break; /* successful completion */

    p = q = r = 0.0;
    if (fabs(e) > tol1) {
      /* fit parabola */
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) p = -p; else q = -q;
      r = e;
      e = d;
    }

    if (fabs(p) < fabs(0.5 * q * r) && p < q * (a - x) && p < q * (b - x)) {
      /* a parabolic interpolation step */
      d = p / q;
      u = x + d;
      /* f must not be evaluated too close to a or b */
      if ((u - a) < t2 || (b - u) < t2)
        d = (x < m) ? tol1 : -tol1;
    } else {
      /* a golden section step */
      e = (x < m) ? b - x : a - x;
      d = GOLDEN * e;
    }

    /* f must not be evaluated too close to x */
    if (fabs(d) >= tol1)
      u = x + d;
    else
      u = x + ((d > 0.0) ? tol1 : -tol1);
    fu = (*f)(u, info);

    /* update a, b, v, w and x */
    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v  = w;  w  = x;  x  = u;
      fv = fw; fw = fx; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v  = w;  w  = u;
        fv = fw; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v  = u;
        fv = fu;
      }
    }
  }

  return x;
}
