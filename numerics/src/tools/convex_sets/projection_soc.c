// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "projection_soc.h"

#include <math.h>    // for sqrt, fabs
#include <stdlib.h>  // for malloc, free

/* Normal part of the projection of x = (xN, xT), t = ||xT||; the tangential part is
   then mu * xN * xT / t on the boundary, 0 in the polar cone. Shared by the 2D, 3D and
   n-dimensional versions (same arithmetic as the former soc_projection_3d). */
static inline unsigned soc_projection_normal(double* xN, double t, double mu) {
  if (mu * t <= -*xN) {
    *xN = 0.0;
    return SOC_PROJECTION_POLAR;
  } else if (t <= mu * *xN) {
    return SOC_PROJECTION_INSIDE;
  } else {
    *xN = (mu * t + *xN) / (mu * mu + 1.0);
    return SOC_PROJECTION_BOUNDARY;
  }
}

unsigned soc_projection_3d(double* x, double mu) {
  /* hypot of libm is safer but really slow */
  double t = sqrt(x[1] * x[1] + x[2] * x[2]);
  unsigned status = soc_projection_normal(&x[0], t, mu);
  if (status == SOC_PROJECTION_BOUNDARY) {
    x[1] = mu * x[0] * x[1] / t;
    x[2] = mu * x[0] * x[2] / t;
  } else if (status == SOC_PROJECTION_POLAR) {
    x[1] = 0.0;
    x[2] = 0.0;
  }
  return status;
}

unsigned soc_projection_2d(double* x, double mu) {
  double t = fabs(x[1]);
  unsigned status = soc_projection_normal(&x[0], t, mu);
  if (status == SOC_PROJECTION_BOUNDARY)
    x[1] = mu * x[0] * ((x[1] > 0.0) ? 1.0 : -1.0);
  else if (status == SOC_PROJECTION_POLAR)
    x[1] = 0.0;
  return status;
}

unsigned soc_projection(double* x, double mu, int dim) {
  if (dim == 3) return soc_projection_3d(x, mu);
  if (dim == 2) return soc_projection_2d(x, mu);
  if (dim == 1) {
    if (x[0] >= 0.0) return SOC_PROJECTION_INSIDE;
    x[0] = 0.0;
    return SOC_PROJECTION_POLAR;
  }
  double t2 = 0.0;
  for (int i = 1; i < dim; i++) t2 += x[i] * x[i];
  double t = sqrt(t2);
  unsigned status = soc_projection_normal(&x[0], t, mu);
  if (status == SOC_PROJECTION_BOUNDARY)
    for (int i = 1; i < dim; i++) x[i] = mu * x[0] * x[i] / t;
  else if (status == SOC_PROJECTION_POLAR)
    for (int i = 1; i < dim; i++) x[i] = 0.0;
  return status;
}

unsigned soc_dual_projection_3d(double* x, double mu) {
  double t = sqrt(x[1] * x[1] + x[2] * x[2]);
  if (t <= -mu * x[0]) {
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    return SOC_PROJECTION_POLAR;
  } else if (mu * t <= x[0]) {
    return SOC_PROJECTION_INSIDE;
  } else {
    double a = (t + mu * x[0]) / (mu * mu + 1.0);
    x[1] = a * x[1] / t;
    x[2] = a * x[2] / t;
    x[0] = mu * a;
    return SOC_PROJECTION_BOUNDARY;
  }
}

unsigned soc_projection_jacobian_3d(double* J, const double* x, double mu) {
  double t = sqrt(x[1] * x[1] + x[2] * x[2]);
  for (int i = 0; i < 9; i++) J[i] = 0.0;

  if (mu * t <= -x[0]) return SOC_PROJECTION_POLAR;
  if (t <= mu * x[0]) {
    J[0] = J[4] = J[8] = 1.0;
    return SOC_PROJECTION_INSIDE;
  }
  /* P = a (1, mu s), a = (xN + mu t)/(1 + mu^2), s = xT/t:
     J = c [1, mu s^T; mu s, mu^2 s s^T] + (mu a/t) [0, 0; 0, I - s s^T], c = 1/(1+mu^2).
     Written as in the former soc_projection_jacobian_3d (same bits). */
  double c = 1. / (1. + mu * mu);
  double s1 = x[1] / t;
  double s2 = x[2] / t;
  double muc = mu * c;
  double muct = muc / t;
  double xNmut = x[0] + mu * t;
  J[0] = c;
  J[1] = muc * s1;
  J[2] = muc * s2;
  J[3] = J[1];
  J[6] = J[2];
  J[4] = muct * (xNmut - x[0] * s1 * s1);
  J[7] = muct * (-x[0] * s1 * s2);
  J[5] = muct * (-x[0] * s1 * s2);
  J[8] = muct * (xNmut - x[0] * s2 * s2);
  return SOC_PROJECTION_BOUNDARY;
}

void soc_projection_3d_stack(int nc, double* x, const double* mu, double* work) {
  int alloc_work = 0;
  if (work == NULL) {
    work = (double*)malloc(nc * sizeof(double));
    alloc_work = 1;
  }
  /* Tangential norms of all the cones at once: a plain loop is faster than per-cone BLAS
     calls on these short contiguous blocks. */
  for (int i = 0; i < nc; ++i) {
    const double t1 = x[3 * i + 1];
    const double t2 = x[3 * i + 2];
    work[i] = sqrt(t1 * t1 + t2 * t2);
  }
  for (int i = 0; i < nc; ++i) {
    /* same arithmetic as the former soc_projection_3d_stack: one scale factor */
    unsigned status = soc_projection_normal(&x[3 * i], work[i], mu[i]);
    if (status == SOC_PROJECTION_BOUNDARY) {
      const double scale = mu[i] * x[3 * i] / work[i];
      x[3 * i + 1] *= scale;
      x[3 * i + 2] *= scale;
    } else if (status == SOC_PROJECTION_POLAR) {
      x[3 * i + 1] = 0.0;
      x[3 * i + 2] = 0.0;
    }
  }
  if (alloc_work) free(work);
}
