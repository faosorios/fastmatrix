/* $ID: linear_iterative.c, last updated 10-14-2021, F.Osorio */

#include "fastmatrix.h"

/* static functions */
static int cg_iter(double *, int, int, double *, double *, int, double);
static int jacobi_iter(double *, int, int, double *, double *, int, double);
static int seidel_iter(double *, int, int, double *, double *, int, double);
/* ..end declarations */

void
cg_solver(double *a, int *lda, int *n, double *b, double *x, int *maxiter, double *tol, int *iter, int *info)
{ /* conjugate gradients method solves the linear system a %*% x = b iteratively,
   * where 'a' is a nonsingular n-by-n matrix, 'b' is the right hand side and 'x'
   * is the approximate solution. */

  /* test the input parameters */
  *info = 0;
  if (*n < 0) {
    *info = -3;
  } else if (*lda < MAX(1, *n)) {
    *info = -2;
  } else if (*maxiter < 0) {
    *info = -6;
  } else if (*tol <= 0.0) {
    *info = -7;
  }
  if (*info != 0) return;

  /* quick return if possible */
  if (*n == 0 || *maxiter == 0)
    return;

  /* call solver */
  *iter = cg_iter(a, *lda, *n, b, x, *maxiter, *tol);
}

static int
cg_iter(double *a, int lda, int n, double *b, double *x, int maxiter, double tol)
{ /* conjugate gradients iteration */
  int iter = 0;
  double alpha, beta, gamma, u, v, *h, *q, *r;
  double accum, az, z, scale, ssq;

  /* initialization */
  h = (double *) Calloc(n, double);
  q = (double *) Calloc(n, double);
  r = (double *) Calloc(n, double);

  /* warming-up */
  scale = 0.0; ssq = 1.0;
  for (int i = 0; i < n; i++) {
    z = q[i] = r[i] = b[i];

    if (z != 0.0) {
      az = fabs(z);

      if (scale < az) {
        ssq = 1.0 + ssq * (scale / az) * (scale / az);
        scale = az;
      } else
        ssq += (az / scale) * (az / scale);
    }
  }
  gamma = SQR(scale * sqrt(ssq));

  /* iteration */
  while (gamma > tol) {
    FM_mult_mat(h, a, lda, n, n, q, n, n, 1);
    /* computing dot product and step-length */
    accum = 0.0; scale = 0.0; ssq = 1.0;
    for (int i = 0; i < n; i++) { /* code re-use! */
      accum += q[i] * h[i];
      z = r[i];

      if (z != 0.0) {
        az = fabs(z);

        if (scale < az) {
          ssq = 1.0 + ssq * (scale / az) * (scale / az);
          scale = az;
        } else
          ssq += (az / scale) * (az / scale);
      }
    }
    u = accum;
    v = SQR(scale * sqrt(ssq));
    alpha = v / u;

    /* approximate solution and residual */
    for (int i = 0; i < n; i++) {
      x[i] += alpha * q[i];
      r[i] -= alpha * h[i];
    }
    /* improvement step */
    beta = FM_norm_sqr(r, 1, n) / v;

    /* updating search direction and convergence criterion (gamma) */
    scale = 0.0; ssq = 1.0;
    for (int i = 0; i < n; i++) {
      z = r[i];
      q[i] = z + beta * q[i];

      if (z != 0.0) {
        az = fabs(z);

        if (scale < az) {
          ssq = 1.0 + ssq * (scale / az) * (scale / az);
          scale = az;
        } else
          ssq += (az / scale) * (az / scale);
      }
    }
    gamma = SQR(scale * sqrt(ssq));

    iter++;
    if (iter > maxiter)
      break; /* maximum number of iterations exceeded */
  }

  Free(h); Free(q); Free(r);

  return(iter);
}

void
jacobi_solver(double *a, int *lda, int *n, double *b, double *x, int *maxiter, double *tol, int *iter, int *info)
{ /* Jacobi's method solves the linear system a %*% x = b iteratively, where 'a'
   * is a nonsingular n-by-n matrix, 'b¡ is the right hand side and 'x' is the
   * approximate solution. */

  /* test the input parameters */
  *info = 0;
  if (*n < 0) {
    *info = -3;
  } else if (*lda < MAX(1, *n)) {
    *info = -2;
  } else if (*maxiter < 0) {
    *info = -6;
  } else if (*tol <= 0.0) {
    *info = -7;
  }
  if (*info != 0) return;

  /* quick return if possible */
  if (*n == 0 || *maxiter == 0)
    return;

  /* check diagonal elements */
  for (int i = 0; i < *n; i++) {
    if (a[i * (*lda + 1)] == 0.0) {
      *info = i + 1;
      return;
    }
  }

  /* call solver */
  *iter = jacobi_iter(a, *lda, *n, b, x, *maxiter, *tol);
}

static int
jacobi_iter(double *a, int lda, int n, double *b, double *x, int maxiter, double tol)
{ /* Jacobi's iteration */
  int iter = 0;
  double accum, az, z, check, scale, ssq, *xnew;

  xnew = (double *) Calloc(n, double);

  /* main loop */
  repeat {
    /* perfom Jacobi iteration */
    for (int i = 0; i < n; i++) {
      accum = 0.0;
      for (int j = 0; j < n; j++) {
        if (i != j)
          accum += a[i + j * lda] * x[j];
      }
      xnew[i] = (b[i] - accum) / a[i * (lda + 1)];
    }

    /* convergence criterion */
    scale = 0.0; ssq = 1.0;
    for (int i = 0; i < n; i++) {
      z = xnew[i] - x[i];
      if (z != 0.0) {
        az = fabs(z);
        if (scale < az) {
          ssq = 1.0 + ssq * (scale / az) * (scale / az);
          scale = az;
        } else
          ssq += (az / scale) * (az / scale);
      }
    }
    check = scale * sqrt(ssq);

    iter++;

    /* eval convergence */
    if (check < tol)
      break; /* successful completion */
    if (iter >= maxiter)
      break; /* maximum number of iterations exceeded */

    /* update solution */
    BLAS1_copy(x, 1, xnew, 1, n);
  }

  Free(xnew);

  return iter;
}

void
seidel_solver(double *a, int *lda, int *n, double *b, double *x, int *maxiter, double *tol, int *iter, int *info)
{ /* Gauss-Seidel method solves the linear system a %*% x = b iteratively, where
   * 'a' is a nonsingular n-by-n matrix, 'b' is the right hand side and 'x' is the
   * approximate solution. */

  /* test the input parameters */
  *info = 0;
  if (*n < 0) {
    *info = -3;
  } else if (*lda < MAX(1, *n)) {
    *info = -2;
  } else if (*maxiter < 0) {
    *info = -6;
  } else if (*tol <= 0.0) {
    *info = -7;
  }
  if (*info != 0) return;

  /* quick return if possible */
  if (*n == 0 || *maxiter == 0)
    return;

  /* check diagonal elements */
  for (int i = 0; i < *n; i++) {
    if (a[i * (*lda + 1)] == 0.0) {
      *info = i + 1;
      return;
    }
  }

  /* call solver */
  *iter = seidel_iter(a, *lda, *n, b, x, *maxiter, *tol);
}

static int
seidel_iter(double *a, int lda, int n, double *b, double *x, int maxiter, double tol)
{ /* Gauss-Seidel iteration */
  int iter = 0;
  double accum, az, z, check, scale, ssq, *xnew;

  xnew = (double *) Calloc(n, double);

  /* main loop */
  repeat {
    /* perfom Gauss-Seidel iteration */
    for (int i = 0; i < n; i++) {
      accum = 0.0;
      for (int j = 0; j < n; j++) {
        if (i < j)
          accum += a[i + j * lda] * x[j];
        else if (i > j)
          accum += a[i + j * lda] * xnew[j];
      }
      xnew[i] = (b[i] - accum) / a[i * (lda + 1)];
    }

    /* convergence criterion */
    scale = 0.0; ssq = 1.0;
    for (int i = 0; i < n; i++) {
      z = xnew[i] - x[i];
      if (z != 0.0) {
        az = fabs(z);
        if (scale < az) {
          ssq = 1.0 + ssq * (scale / az) * (scale / az);
          scale = az;
        } else
          ssq += (az / scale) * (az / scale);
      }
    }
    check = scale * sqrt(ssq);

    iter++;

    /* eval convergence */
    if (check < tol)
      break; /* successful completion */
    if (iter >= maxiter)
      break; /* maximum number of iterations exceeded */

    /* update solution */
    BLAS1_copy(x, 1, xnew, 1, n);
  }

  Free(xnew);

  return iter;
}
