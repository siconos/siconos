// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "natural_map_error.h"

#include <math.h>    // for sqrt, fabs
#include <stddef.h>  // for NULL

#include "SolverOptions.h"  // for SolverOptions

static inline double safeguard(double s) { return (s > 1.0) ? s : 1.0; }

/* the norms are accumulated block by block: sums of squares for the Euclidean and
   weighted norms, sums for the cone norms */
static inline int norm_is_quadratic(int norm) {
  return norm == NATURAL_MAP_ERROR_NORM_EUCLIDEAN ||
         norm == NATURAL_MAP_ERROR_NORM_CONE_WEIGHTED;
}

/* contribution of one block x = (x_N, x_T) to the accumulated norm */
static inline double block_norm_term(int norm, const double* x, int dim, double mu,
                                     double alpha) {
  double t2 = 0.0;
  for (int k = 1; k < dim; k++) t2 += x[k] * x[k];
  switch (norm) {
    case NATURAL_MAP_ERROR_NORM_EUCLIDEAN:
      return x[0] * x[0] + t2;
    case NATURAL_MAP_ERROR_NORM_CONE:
      return fabs(x[0]) + sqrt(t2);
    case NATURAL_MAP_ERROR_NORM_CONE_MU:
      return fabs(x[0]) + mu * sqrt(t2);
    case NATURAL_MAP_ERROR_NORM_CONE_WEIGHTED:
      return x[0] * x[0] + ((alpha > 0.0) ? alpha : 0.0) * t2;
    default:
      return 0.0;
  }
}

static inline double finish_norm(int norm, double acc) {
  return norm_is_quadratic(norm) ? sqrt(acc) : acc;
}

static inline double reduce(int reduction, double nr, double nv) {
  switch (reduction) {
    case NATURAL_MAP_ERROR_REDUCTION_MEAN:
      return sqrt(safeguard(nr) * safeguard(nv));
    case NATURAL_MAP_ERROR_REDUCTION_MIXED:
      return safeguard(nr + nv);
    case NATURAL_MAP_ERROR_REDUCTION_MAX:
    default: {
      const double sr = safeguard(nr), sv = safeguard(nv);
      return (sr < sv) ? sr : sv;
    }
  }
}

/* coefficient of a block for the CONE_MU norm: mu for the cones, 1 otherwise */
static inline double cone_coefficient(const ConvexSet* K, int i) {
  return ((K->type == CONVEX_SET_SOC || K->type == CONVEX_SET_ROLLING_CONE) && K->p1)
             ? K->p1[i]
             : 1.0;
}

NaturalMapErrorScaling natural_map_error_scaling(const SolverOptions* options,
                                                 NaturalMapErrorScaling defaults) {
  if (!options) return defaults;
  NaturalMapErrorScaling s;
  s.norm = options->iparam[NATURAL_MAP_ERROR_IPARAM_NORM];
  s.scope = options->iparam[NATURAL_MAP_ERROR_IPARAM_SCOPE];
  s.reduction = options->iparam[NATURAL_MAP_ERROR_IPARAM_REDUCTION];
  s.alpha = options->dparam[NATURAL_MAP_ERROR_DPARAM_ALPHA];
  return s;
}

double natural_map_error_norm(int norm, const double* x, int m, int dim, const double* mu,
                              double alpha) {
  if (norm == NATURAL_MAP_ERROR_NORM_NO) return 1.0;
  if (m <= 0 || dim <= 0 || x == NULL) return 0.0;
  double acc = 0.0;
  for (int i = 0; i < m; i++)
    acc += block_norm_term(norm, x + i * dim, dim, mu ? mu[i] : 1.0, alpha);
  return finish_norm(norm, acc);
}

double natural_map_error_scale(int norm, int reduction, const double* x, const double* y,
                               int m, int dim, const double* mu, double alpha) {
  if (norm == NATURAL_MAP_ERROR_NORM_NO) return 1.0;
  return reduce(reduction, natural_map_error_norm(norm, x, m, dim, mu, alpha),
                natural_map_error_norm(norm, y, m, dim, mu, alpha));
}

void natural_map_error_init(NaturalMapErrorSum* sum) {
  sum->F2 = 0.0;
  sum->nr = 0.0;
  sum->nv = 0.0;
}

/* ||r - P_K(r - v)||^2 for a cone of R^3, inlined (hot path of fc3d, gfc3d, Drucker-Prager):
   same arithmetic as soc_projection_3d. */
static inline double soc3_natural_map_F2(const double* r, const double* v, double mu) {
  const double w0 = r[0] - v[0], w1 = r[1] - v[1], w2 = r[2] - v[2];
  const double t = sqrt(w1 * w1 + w2 * w2);
  double p0, p1, p2;
  if (mu * t <= -w0) {
    p0 = p1 = p2 = 0.0;
  } else if (t <= mu * w0) {
    p0 = w0;
    p1 = w1;
    p2 = w2;
  } else {
    p0 = (mu * t + w0) / (mu * mu + 1.0);
    p1 = mu * p0 * w1 / t;
    p2 = mu * p0 * w2 / t;
  }
  const double f0 = r[0] - p0, f1 = r[1] - p1, f2 = r[2] - p2;
  return f0 * f0 + f1 * f1 + f2 * f2;
}

/* Loop over the blocks; F2_BLOCK(i, ri, vi) gives ||F_i||^2. Written as a macro so that the
   SOC 3D loop is compiled with the inlined projection and the generic one with the call. */
#define NATURAL_MAP_ERROR_LOOP(F2_BLOCK, DIM)                                                 \
  for (int i = 0; i < K->nblocks; i++) {                                                      \
    const double* ri = r + i * (DIM);                                                         \
    const double* vi = v + i * (DIM);                                                         \
    double F2 = F2_BLOCK(i, ri, vi);                                                          \
    if (scaled) {                                                                             \
      const double mu = cone_coefficient(K, i);                                               \
      const double tr = block_norm_term(norm, ri, (DIM), mu, s->alpha);                       \
      const double tv = block_norm_term(norm, vi, (DIM), mu, s->alpha);                       \
      if (per_block) {                                                                        \
        const double sc = reduce(s->reduction, finish_norm(norm, tr), finish_norm(norm, tv)); \
        F2 /= sc * sc;                                                                        \
      } else {                                                                                \
        sum->nr += tr;                                                                        \
        sum->nv += tv;                                                                        \
      }                                                                                       \
    }                                                                                         \
    sum->F2 += F2;                                                                            \
  }

#define SOC3_F2(i, ri, vi) soc3_natural_map_F2(ri, vi, K->p1[i])
#define GENERIC_F2(i, ri, vi) convex_set_natural_map_block(K, i, ri, vi, NULL)

void natural_map_error_add(const NaturalMapErrorScaling* s, const ConvexSet* K,
                           const double* r, const double* v, NaturalMapErrorSum* sum) {
  const int dim = K->dim;
  const int norm = s->norm;
  const int scaled = (norm != NATURAL_MAP_ERROR_NORM_NO);
  const int per_block = scaled && (s->scope == NATURAL_MAP_ERROR_SCOPE_PER_BLOCK);

  if (K->type == CONVEX_SET_SOC && dim == 3) {
    /* literal dimension: the norm loops are unrolled */
    NATURAL_MAP_ERROR_LOOP(SOC3_F2, 3)
  } else {
    NATURAL_MAP_ERROR_LOOP(GENERIC_F2, dim)
  }
}

#undef SOC3_F2
#undef GENERIC_F2
#undef NATURAL_MAP_ERROR_LOOP

double natural_map_error_value(const NaturalMapErrorScaling* s,
                               const NaturalMapErrorSum* sum) {
  double e = sqrt(sum->F2);
  if (s->norm != NATURAL_MAP_ERROR_NORM_NO && s->scope == NATURAL_MAP_ERROR_SCOPE_PER_VECTOR)
    e /= reduce(s->reduction, finish_norm(s->norm, sum->nr), finish_norm(s->norm, sum->nv));
  return e;
}

double natural_map_error(const NaturalMapErrorScaling* s, const ConvexSet* K, const double* r,
                         const double* v) {
  NaturalMapErrorSum sum;
  natural_map_error_init(&sum);
  natural_map_error_add(s, K, r, v, &sum);
  return natural_map_error_value(s, &sum);
}
