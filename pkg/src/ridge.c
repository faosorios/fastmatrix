/* $ID: ridge.c, last updated 2020-12-27, F.Osorio */

#include "fastmatrix.h"
#include "ridge.h"

/* static functions */
static void ridge_default(double, GCVinfo);
static void ridge_grid(double *, int, double *, GCVinfo, double *);
static void ridge_GCV(double *, GCVinfo, double);
static double log_GCV(double, void *);
/* ..end declarations */

void
OLS_ridge(double *x, int *ldx, int *nrow, int *ncol, double *y, double *coef, double *scale,
  double *fitted, double *resid, double *rss, double *edf, double *pen, double *gcv,
  double *hkb, double *lw, double *lambda, double *lambda_opt, int *ngrid, int *task,
  double *tolerance)
{ /* ridge regression */
  int errcode = 0, job, n = *nrow, p = *ncol;
  double *a, *d, *v, *rhs, s2;
  double HKB, LW, PEN;
  GCVinfo pars;

  a    = (double *) Calloc(p, double);
  d    = (double *) Calloc(p, double);
  rhs  = (double *) Calloc(p, double);
  v    = (double *) Calloc(p * p, double);
  pars = (GCVinfo)  Calloc(1, GCV_info);

  /* SVD of the model matrix */
  job = 21; /* left singular vectors overwite 'x' */
  FM_svd_decomp(x, *ldx, n, p, x, *ldx, d, v, p, job, &errcode);
  if (errcode)
    error("Singular value decomposition gave code %d", errcode);

  /* compute the right-hand side */
  FM_crossprod(rhs, x, *ldx, n, p, y, n, n, 1);

  /* HKB and LW estimates */
  for (int j = 0; j < p; j++)
    a[j]  = rhs[j] / d[j];
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

  Free(a); Free(d); Free(rhs); Free(v); Free(pars);
}

static void
ridge_default(double lambda, GCVinfo st)
{
  int n = st->n, p = st->p;
  double *a, div, edf = 0.0, s2;

  a = (double *) Calloc(p, double);

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

  Free(a);
}

static void
ridge_grid(double *lambda, int ngrid, double *gcv, GCVinfo st, double *lambda_opt)
{
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
{
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

  a = (double *) Calloc(p, double);

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

  Free(a);

  return val;
}
