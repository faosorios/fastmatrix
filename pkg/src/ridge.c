/* $ID: ridge.c, last updated 2024-09-03, F.Osorio */

#include "fastmatrix.h"
#include "ridge.h"

/* static functions */
static void ridge_default(double, GCVinfo);
static void ridge_grid(double *, int, double *, GCVinfo, double *);
static void ridge_GCV(double *, GCVinfo, double);
static double log_GCV(double, void *);
static double fnc1_dot(double, double, double, double);
static double fnc1_ddot(double, double, double, double);
static void ridge_ORP1(double *, GCVinfo, double, int *);
/* ..end declarations */

void
OLS_ridge(double *x, int *ldx, int *nrow, int *ncol, double *y, double *coef, double *scale,
  double *fitted, double *resid, double *rss, double *edf, double *pen, double *gcv,
  double *hkb, double *lw, double *lambda, double *lambda_opt, int *ngrid, int *task,
  double *tolerance, int *maxiter)
{ /* ridge regression */
  int errcode = 0, job, n = *nrow, p = *ncol;
  double *a, *d, *v, *rhs, s2;
  double HKB, LW, PEN;
  GCVinfo pars;

  a    = (double *) R_Calloc(p, double);
  d    = (double *) R_Calloc(p, double);
  rhs  = (double *) R_Calloc(p, double);
  v    = (double *) R_Calloc(p * p, double);
  pars = (GCVinfo)  R_Calloc(1, GCV_info);

  /* SVD of the model matrix */
  job = 21; /* left singular vectors overwite 'x' */
  FM_svd_decomp(x, *ldx, n, p, x, *ldx, d, v, p, job, &errcode);
  if (errcode)
    error("Singular value decomposition gave code %d", errcode);

  /* compute the right-hand side */
  FM_crossprod(rhs, x, *ldx, n, p, y, n, n, 1);

  /* HKB and LW estimates */
  for (int j = 0; j < p; j++)
    a[j] = rhs[j] / d[j];
  PEN = FM_norm_sqr(a, 1, p);
  for (int j = 0; j < p; j++)
    a[j] *= d[j];
  FM_mult_mat(fitted, x, *ldx, n, p, a, p, p, 1);
  for (int i = 0; i < n; i++)
    resid[i] = y[i] - fitted[i];
  s2 = FM_norm_sqr(resid, 1, n) / (n - p);
  HKB = p * s2 / PEN;
  LW  = n * p * s2 / FM_norm_sqr(fitted, 1, n);

  /* constructs a GCV object */
  pars->n = n;
  pars->p = p;
  pars->u = x;
  pars->d = d;
  pars->y = y;
  pars->rhs = rhs;
  pars->a = a;
  pars->fitted = fitted;
  pars->resid = resid;

  switch (*task) {
    case 0: /* none */
      ridge_default(*lambda, pars);
      break;
    case 1: /* grid */
      ridge_grid(lambda, *ngrid, gcv, pars, lambda_opt);
      ridge_default(*lambda_opt, pars); /* running with 'optimal' lambda */
      break;
    case 2: /* GCV */
      ridge_GCV(lambda, pars, *tolerance);
      break;
    case 3: /* ORP1 */
      ridge_ORP1(lambda, pars, *tolerance, maxiter);
      ridge_default(*lambda, pars); /* running with 'optimal' lambda */
      break;
    default:
      ridge_default(*lambda, pars);
      break;
  }
  FM_mult_mat(coef, v, p, p, p, pars->a, p, p, 1);

  /* saving results */
  *scale = (pars->RSS + *lambda * pars->pen) / n;
  *edf = pars->edf;
  *pen = pars->pen;
  *gcv = pars->GCV; /* actually log(GCV) */
  *hkb = HKB;
  *lw  = LW;
  *rss = pars->RSS;

  R_Free(a); R_Free(d); R_Free(rhs); R_Free(v); R_Free(pars);
}

static void
ridge_default(double lambda, GCVinfo st)
{ /* computes ridge estimator (using fixed regulation parameter) */
  int n = st->n, p = st->p;
  double *a, div, edf = 0.0, s2;

  a = (double *) R_Calloc(p, double);

  /* compute the coefficients and effective degrees of freedom */
  for (int j = 0; j < p; j++) {
    div  = lambda + SQR((st->d)[j]);
    edf += SQR((st->d)[j]) / div;
    (st->a)[j] = (st->d)[j] * (st->rhs)[j] / div;
  }
  st->pen = FM_norm_sqr(st->a, 1, p);

  /* compute fitted values and residuals */
  for (int j = 0; j < p; j++)
    a[j] = (st->a)[j] * (st->d)[j];
  FM_mult_mat(st->fitted, st->u, n, n, p, a, p, p, 1);
  for (int i = 0; i < n; i++)
    (st->resid)[i] = (st->y)[i] - (st->fitted)[i];

  /* compute GCV criterion */
  st->RSS = FM_norm_sqr(st->resid, 1, n);
  s2 = st->RSS / (n - edf);
  st->edf = edf;
  st->GCV = s2 / (1.0 - edf / n);

  R_Free(a);
}

static void
ridge_grid(double *lambda, int ngrid, double *gcv, GCVinfo st, double *lambda_opt)
{ /* select optimal ridge (regulation) parameter over a grid */
  double GCV_min, opt = 0.0;

  /* compute GCV criterion over a grid */
  ridge_default(lambda[0], st);
  gcv[0] = GCV_min = st->GCV;
  for (int k = 1; k < ngrid; k++) {
    ridge_default(lambda[k], st);
    gcv[k] = st->GCV;
    if (gcv[k] <= GCV_min) {
      GCV_min = gcv[k];
      opt = lambda[k];
    }
  }
  *lambda_opt = opt;
}

static void
ridge_GCV(double *lambda, GCVinfo st, double tol)
{ /* select optimal ridge parameter minimizing the GCV criterion */
  double conv, upper_lambda;
  const double phi = 2.0 - GOLDEN;

  /* call optimizer */
  upper_lambda = *lambda;
  do {
    *lambda = brent(0., upper_lambda, log_GCV, st, tol);
    conv = fabs(*lambda - upper_lambda);
    upper_lambda *= phi;
  } while (conv < tol);
}

double
log_GCV(double lambda, void *pars)
{ /* for brent's procedure */
  GCVinfo st = (GCVinfo) pars;
  int n = st->n, p = st->p;
  double *a, div, edf = 0.0, s2, val;

  a = (double *) R_Calloc(p, double);

  /* compute the coefficients and effective degrees of freedom */
  for (int j = 0; j < p; j++) {
    div  = lambda + SQR((st->d)[j]);
    edf += SQR((st->d)[j]) / div;
    (st->a)[j] = (st->d)[j] * (st->rhs)[j] / div;
  }
  st->pen = FM_norm_sqr(st->a, 1, p);

  /* compute fitted values and residuals */
  for (int j = 0; j < p; j++)
    a[j] = (st->a)[j] * (st->d)[j];
  FM_mult_mat(st->fitted, st->u, n, n, p, a, p, p, 1);
  for (int i = 0; i < n; i++)
    (st->resid)[i] = (st->y)[i] - (st->fitted)[i];

  /* compute log-GCV criterion */
  st->RSS = FM_norm_sqr(st->resid, 1, n);
  s2  = st->RSS / (n - edf);
  val = log(s2) - log(1.0 - edf / n);
  st->edf = edf;
  st->GCV = val;

  R_Free(a);

  return val;
}

/* functions to evaluate the MSE criterion and its derivatives, which are called
 * by ridge_ORP1 (nested functions are forbidden in ISO C) */
static double
fnc1_dot(double alpha, double d2, double s2, double k) {
  return d2 * (k * SQR(alpha) - s2) / CUBE(d2 + k);
}
static double
fnc1_ddot(double alpha, double d2, double s2, double k) {
  return d2 * (SQR(alpha) * (d2 - 2.0 * k) + 3.0 * s2) / R_pow_di(d2 + k, 4);
}

static void
ridge_ORP1(double *lambda, GCVinfo st, double tol, int *maxit)
{ /* select optimal ridge parameter minimizing the mean squared estimation (MSE) error 
   * based on AS 223: Applied Statistics 36, 1987, 112-118. doi: 10.2307/2347851 */
  int n = st->n, p = st->p, iter = 0;
  double check, d2, f1_dot, f1_ddot, k, knew, s2;

  /* compute s2 estimator */
  s2 = FM_norm_sqr(st->resid, 1, n) / (n - p);

  /* initialization */
  k = *lambda;

  /* main loop */
  repeat {
    f1_dot = f1_ddot = 0.0;
    for (int j = 0; j < p; j++) {
      d2 = SQR((st->d)[j]);
      f1_dot += fnc1_dot((st->a)[j], d2, s2, k);
      f1_ddot += fnc1_ddot((st->a)[j], d2, s2, k);
    }

    knew = k - f1_dot / f1_ddot;
    check = fabs(knew - k);
    iter++;

    /* eval convergence */
    if (check < tol)
      break; /* successful completion */
    if (iter >= *maxit)
      break; /* maximum number of iterations exceeded */

    /* update solution */
    k = knew;
  }

  *lambda = knew;
  *maxit = iter;
}
