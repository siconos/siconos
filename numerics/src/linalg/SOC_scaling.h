/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SOC_scaling_H
#define SOC_scaling_H

/*!\file SOC_scaling.h
  \brief Norms and scale functions for stacked second-order cone problems.

  The functions operate on a vector \f$ x \in \mathbb{R}^{m \cdot dim} \f$
  made of \p m cones of dimension \p dim.  Single-vector norms
  (`SOC_norm_*`) are provided, together with three two-vector scalings:
  - mixed scalings (`SOC_mixed_scaling_*`) obtained by adding the norms of
    two vectors;
  - per-variable scalings (`SOC_scaling_*`) obtained by taking the maximum of
    the two norms;
  - geometric per-variable scalings (`SOC_geom_scaling_*`) obtained by taking
    the geometric mean of the two norms.

  Safeguarding (`max(1, sigma)`) is intentionally left to the caller so that
  it is applied in one central place.
*/

#if defined(__cplusplus)
extern "C" {
#endif

/* -------------------------------------------------------------------------- */
/* Single-vector norms                                                        */
/* -------------------------------------------------------------------------- */

/**
    Euclidean norm of a stacked cone vector.

    \f$ \|x\|_2 = \sqrt{\sum_{k=0}^{m \cdot dim -1} x_k^2} \f$

    \param[in] x   vector of length m * dim
    \param[in] m   number of cones
    \param[in] dim dimension of one cone
    \return the Euclidean norm
*/
double SOC_norm_euclidean(const double* x, int m, int dim);

/**
    Cone (spectral/Jordan) norm of a stacked cone vector.

    For each cone \f$ j \f$:
    \f$ \|x_j\|_K = |x_{j \cdot dim}| + \sqrt{\sum_{i=1}^{dim-1} x_{j \cdot dim+i}^2} \f$
    and the returned value is \f$ \sum_{j=0}^{m-1} \|x_j\|_K \f$.

    \param[in] x   vector of length m * dim
    \param[in] m   number of cones
    \param[in] dim dimension of one cone
    \return the cone norm
*/
double SOC_norm_cone(const double* x, int m, int dim);

/**
    Friction-cone norm of a stacked cone vector.

    For each cone \f$ j \f$ with friction coefficient \f$ \mu_j \f$:
    \f$ \|x_j\|_{K_{\mu_j}} = |x_{j \cdot dim}| +
        \mu_j \sqrt{\sum_{i=1}^{dim-1} x_{j \cdot dim+i}^2} \f$
    and the returned value is \f$ \sum_{j=0}^{m-1} \|x_j\|_{K_{\mu_j}} \f$.

    \param[in] x   vector of length m * dim
    \param[in] m   number of cones
    \param[in] mu  friction coefficients, array of length m
    \param[in] dim dimension of one cone
    \return the friction-cone norm
*/
double SOC_norm_cone_mu(const double* x, int m, const double* mu, int dim);

/**
    Weighted cone norm of a stacked cone vector.

    \f$ \|x\|_{K,\alpha} = \sqrt{\sum_{j=0}^{m-1}
        \left( x_{j \cdot dim}^2 + \alpha \sum_{i=1}^{dim-1} x_{j \cdot dim+i}^2
        \right)} \f$

    \param[in] x     vector of length m * dim
    \param[in] m     number of cones
    \param[in] dim   dimension of one cone
    \param[in] alpha weight on the tangential part (must be positive)
    \return the weighted cone norm
*/
double SOC_norm_cone_weighted(const double* x, int m, int dim, double alpha);

/* -------------------------------------------------------------------------- */
/* Mixed scalings (sum of the corresponding single-vector norms)              */
/* -------------------------------------------------------------------------- */

/**
    Euclidean mixed scale.

    \f$ \sigma^\mathrm{mixed}_E(x,y) = \|x\|_2 + \|y\|_2 \f$

    \param[in] x   first vector of length m * dim
    \param[in] y   second vector of length m * dim
    \param[in] m   number of cones
    \param[in] dim dimension of one cone
    \return the Euclidean mixed scale
*/
double SOC_mixed_scaling_euclidean(const double* x, const double* y, int m, int dim);

/**
    Cone (spectral/Jordan) mixed scale.

    \f$ \sigma^\mathrm{mixed}_K(x,y) = \|x\|_K + \|y\|_K \f$

    \param[in] x   first vector of length m * dim
    \param[in] y   second vector of length m * dim
    \param[in] m   number of cones
    \param[in] dim dimension of one cone
    \return the cone mixed scale
*/
double SOC_mixed_scaling_cone(const double* x, const double* y, int m, int dim);

/**
    Friction-cone mixed scale.

    \f$ \sigma^\mathrm{mixed}_{K_\mu}(x,y) = \|x\|_{K_\mu} + \|y\|_{K_\mu} \f$

    \param[in] x   first vector of length m * dim
    \param[in] y   second vector of length m * dim
    \param[in] m   number of cones
    \param[in] mu  friction coefficients, array of length m
    \param[in] dim dimension of one cone
    \return the friction-cone mixed scale
*/
double SOC_mixed_scaling_cone_mu(const double* x, const double* y, int m, const double* mu,
                                 int dim);

/**
    Weighted cone mixed scale.

    \f$ \sigma^\mathrm{mixed}_{K,\alpha}(x,y) =
        \|x\|_{K,\alpha} + \|y\|_{K,\alpha} \f$

    \param[in] x     first vector of length m * dim
    \param[in] y     second vector of length m * dim
    \param[in] m     number of cones
    \param[in] dim   dimension of one cone
    \param[in] alpha weight on the tangential part (must be positive)
    \return the weighted cone mixed scale
*/
double SOC_mixed_scaling_cone_weighted(const double* x, const double* y, int m, int dim,
                                       double alpha);

/* -------------------------------------------------------------------------- */
/* Per-variable scalings (max of the corresponding single-vector norms)       */
/* -------------------------------------------------------------------------- */

/**
    Euclidean per-variable scale.

    \f$ \sigma_E(x,y) = \max(\|x\|_2, \|y\|_2) \f$

    \param[in] x   first vector of length m * dim
    \param[in] y   second vector of length m * dim
    \param[in] m   number of cones
    \param[in] dim dimension of one cone
    \return the Euclidean per-variable scale
*/
double SOC_scaling_euclidean(const double* x, const double* y, int m, int dim);

/**
    Cone (spectral/Jordan) per-variable scale.

    \f$ \sigma_K(x,y) = \max(\|x\|_K, \|y\|_K) \f$

    \param[in] x   first vector of length m * dim
    \param[in] y   second vector of length m * dim
    \param[in] m   number of cones
    \param[in] dim dimension of one cone
    \return the cone per-variable scale
*/
double SOC_scaling_cone(const double* x, const double* y, int m, int dim);

/**
    Friction-cone per-variable scale.

    \f$ \sigma_{K_\mu}(x,y) = \max(\|x\|_{K_\mu}, \|y\|_{K_\mu}) \f$

    \param[in] x   first vector of length m * dim
    \param[in] y   second vector of length m * dim
    \param[in] m   number of cones
    \param[in] mu  friction coefficients, array of length m
    \param[in] dim dimension of one cone
    \return the friction-cone per-variable scale
*/
double SOC_scaling_cone_mu(const double* x, const double* y, int m, const double* mu, int dim);

/**
    Weighted cone per-variable scale.

    \f$ \sigma_{K,\alpha}(x,y) =
        \max(\|x\|_{K,\alpha}, \|y\|_{K,\alpha}) \f$

    \param[in] x     first vector of length m * dim
    \param[in] y     second vector of length m * dim
    \param[in] m     number of cones
    \param[in] dim   dimension of one cone
    \param[in] alpha weight on the tangential part (must be positive)
    \return the weighted cone per-variable scale
*/
double SOC_scaling_cone_weighted(const double* x, const double* y, int m, int dim,
                                 double alpha);

/* -------------------------------------------------------------------------- */
/* Geometric per-variable scalings                                            */
/* -------------------------------------------------------------------------- */

/**
    Euclidean geometric per-variable scale.

    \f$ \sigma^\mathrm{geom}_E(x,y) = \sqrt{\|x\|_2 \|y\|_2} \f$

    \param[in] x   first vector of length m * dim
    \param[in] y   second vector of length m * dim
    \param[in] m   number of cones
    \param[in] dim dimension of one cone
    \return the Euclidean geometric per-variable scale
*/
double SOC_geom_scaling_euclidean(const double* x, const double* y, int m, int dim);

/**
    Cone (spectral/Jordan) geometric per-variable scale.

    \f$ \sigma^\mathrm{geom}_K(x,y) = \sqrt{\|x\|_K \|y\|_K} \f$

    \param[in] x   first vector of length m * dim
    \param[in] y   second vector of length m * dim
    \param[in] m   number of cones
    \param[in] dim dimension of one cone
    \return the cone geometric per-variable scale
*/
double SOC_geom_scaling_cone(const double* x, const double* y, int m, int dim);

/**
    Friction-cone geometric per-variable scale.

    \f$ \sigma^\mathrm{geom}_{K_\mu}(x,y) =
        \sqrt{\|x\|_{K_\mu} \|y\|_{K_\mu}} \f$

    \param[in] x   first vector of length m * dim
    \param[in] y   second vector of length m * dim
    \param[in] m   number of cones
    \param[in] mu  friction coefficients, array of length m
    \param[in] dim dimension of one cone
    \return the friction-cone geometric per-variable scale
*/
double SOC_geom_scaling_cone_mu(const double* x, const double* y, int m, const double* mu,
                                int dim);

/**
    Weighted cone geometric per-variable scale.

    \f$ \sigma^\mathrm{geom}_{K,\alpha}(x,y) =
        \sqrt{\|x\|_{K,\alpha} \|y\|_{K,\alpha}} \f$

    \param[in] x     first vector of length m * dim
    \param[in] y     second vector of length m * dim
    \param[in] m     number of cones
    \param[in] dim   dimension of one cone
    \param[in] alpha weight on the tangential part (must be positive)
    \return the weighted cone geometric per-variable scale
*/
double SOC_geom_scaling_cone_weighted(const double* x, const double* y, int m, int dim,
                                      double alpha);

#if defined(__cplusplus)
}
#endif

#endif
