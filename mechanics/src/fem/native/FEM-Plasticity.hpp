// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! @file FEM-Plasticity.hpp
  One-step plasticity problem for SolidLinearTIDS integrated with
  fem::integrators::MoreauJeanGOSI, solved by condensation on the constrained stresses and
  the numerics PlasticityProblem.
*/
#ifndef FEMPlasticity_H
#define FEMPlasticity_H

#include "FEM-GlobalFrictionContact.hpp"
#include "Plasticity_options.h"  // for PLASTICITY_2D_NSGS

namespace siconos::mechanics::fem::nonsmooth_formulations {

/** One-step plasticity problem of solids with stress unknowns.

    The global problem assembled by GlobalFrictionContact (iteration matrix W of
    fem::integrators::MoreauJeanGOSI, boundary conditions included) reads

    \f[ W x = q + H z, \quad y = H^T x + b, \quad 0 \in z + N_K(y) \f]

    with \f$ x = [v; \sigma] \f$, y the stresses seen by the plasticity interactions
    (fem::StressLinearTIR, \f$ y = C\sigma \f$), K the elastic domain of their law and
    \f$ z = \lambda \f$ (\f$ C^T\lambda = -h\dot\epsilon_p \f$, associated flow rule).

    Let c be the rows of x touched by H (the constrained stresses) and r the others
    (velocities, unconstrained stresses). H restricted to c, \f$ H_c \f$, must be square
    and invertible (each constrained stress belongs to exactly one interaction, with an
    invertible C block). Eliminating \f$ x_r \f$,

    \f[ W_{red} = W_{cc} - W_{cr} W_{rr}^{-1} W_{rc},\quad
        \tilde q = q_c - W_{cr} W_{rr}^{-1} q_r,\quad W_{red}\, x_c = \tilde q + H_c z, \f]

    and with \f$ x_c = H_c^{-T}(y - b) \f$:

    \f[ z = M_p\, y + q_p,\quad M_p = H_c^{-1} W_{red} H_c^{-T},\quad
        q_p = -M_p b - H_c^{-1}\tilde q, \f]

    which is the numerics PlasticityProblem \f$ 0 \in M_p y + q_p + N_K(y) \f$ (y plays
    the role of the stress, z of the strain rate). When \f$ W_{rr} \f$ is diagonal (lumped
    mass and all stresses constrained) it is inverted entry-wise and \f$ W_{red} \f$ is
    sparse; otherwise a sparse LU of \f$ W_{rr} \f$ is used.

    Supported law: VonMisesPlasticityNSL (Voigt coordinates, the hypothesis of the law;
    all the interactions must share it). y is the physical stress when C is the selection
    of the element stresses divided by sqrt(area * thickness) (see
    SolidLinearTIDS::elementStressSelection).

    See docs/devel/FEM_PLASTICITY_SCHEME.md.
*/
class Plasticity : public GlobalFrictionContact {
 protected:
  ACCEPT_SERIALIZATION(Plasticity);

  /** solve the condensed plasticity problem, fill _globalVelocities, _z (multipliers)
   *  and _w (y)
   *  \return the numerics solver info
   */
  int solve_condensed();

 public:
  /** constructor
   *  \param numericsSolverId id of the numerics plasticity solver
   */
  Plasticity(int numericsSolverId = PLASTICITY_2D_NSGS);

  /** destructor */
  virtual ~Plasticity() noexcept = default;

  /** compute the one-step problem (assembly by GlobalFrictionContact::preCompute,
   *  condensation, numerics solve, postCompute)
   *  \param time current time
   *  \return the numerics solver info
   */
  int compute(double time) override;
};

}  // namespace siconos::mechanics::fem::nonsmooth_formulations

#endif
