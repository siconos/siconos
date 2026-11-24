/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "SphereNEDS.hpp"

#include "SiconosVector.hpp"
#include "StorageTools.hpp"

siconos::collision::native::bodies::SphereNEDS::SphereNEDS(
    double r, double m, Eigen::Ref<siconos::algebra::SiconosMatrix33> inertia,
    Eigen::Ref<siconos::algebra::SiconosVector7> qinit,
    Eigen::Ref<siconos::algebra::SiconosVector6> vinit, siconos::algebra::AliasTag)
    : siconos::modeling::NewtonEulerDS{qinit, vinit, m, inertia, siconos::algebra::alias_t},
      radius{r} {}

siconos::collision::native::bodies::SphereNEDS::SphereNEDS(
    double r, double m, const siconos::algebra::SiconosMatrix33& inertia,
    const siconos::algebra::SiconosVector7& qinit,
    const siconos::algebra::SiconosVector6& vinit, siconos::algebra::CopyTag)
    : siconos::modeling::NewtonEulerDS{qinit, vinit, m, inertia, siconos::algebra::copy_t},
      radius{r} {}

double siconos::collision::native::bodies::SphereNEDS::getQ(siconos::algebra::Index pos) {
  assert(pos < 7);
  return ((*state_q_)(pos));
};

double siconos::collision::native::bodies::SphereNEDS::getTwist(siconos::algebra::Index pos) {
  assert(pos < 6);
  return ((*twist_)(pos));
};
