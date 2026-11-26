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

#include <string>
#include <vector>

#include "FEM-GlobalFrictionContact.hpp"
#include "Plasticity_options.h"  // for PLASTICITY_2D_NSN

namespace siconos::mechanics::fem::nonsmooth_formulations {

/** condensed problem, built once and reused (FEM-Plasticity.cpp) */
struct PlasticityCondensation;

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
    the role of the stress, z of the strain rate), solved by default by the global
    semi-smooth Newton PLASTICITY_2D_NSN. When \f$ W_{rr} \f$ is diagonal (lumped
    mass and all stresses constrained) it is inverted entry-wise and \f$ W_{red} \f$ is
    sparse; otherwise a sparse LU of \f$ W_{rr} \f$ is used.

    Supported laws (one model for all the interactions; numerics Voigt coordinates):
    VonMisesPlasticityNSL (the hypothesis of the law, shared by all the interactions) and
    DruckerPragerPlasticityNSL (cone slope, dilatancy and cohesion per interaction). y is the
   physical stress when C is the selection of the element stresses divided by sqrt(area *
   thickness) (see SolidLinearTIDS::elementStressSelection).

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

  /** build the condensed problem from the assembled W, H (GlobalFrictionContact::preCompute):
   *  selections, Schur complement, M_p, numerics problem (model parameters) */
  void assemble();

  /** the condensed problem, built once (the iteration matrix and the relations are
   *  time-invariant); rebuilt when the number of interactions changes */
  std::shared_ptr<PlasticityCondensation> _condensation{nullptr};

  /** number of calls of compute with a non-empty problem (one per time step) */
  unsigned int _n_compute{0};

  /** prefix of the files of the dumped numerics problems (no dump if empty) */
  std::string _dump_prefix{};

  /** calls of compute (1, 2, ...) whose problem is dumped (all if empty) */
  std::vector<unsigned int> _dump_steps{};

 public:
  /** constructor
   *  \param numericsSolverId id of the numerics plasticity solver: PLASTICITY_2D_NSN
   *  (global semi-smooth Newton, default: tens of iterations whatever the conditioning of
   *  the condensed matrix) or PLASTICITY_2D_NSGS (cost O(cond) sweeps in quasi-static
   *  regimes)
   */
  Plasticity(int numericsSolverId = PLASTICITY_2D_NSN);

  /** destructor */
  virtual ~Plasticity() noexcept;

  /** assembly: the full problem (W, H, q, b) and its condensation at the first call (or when
   *  the interactions change), then only the right-hand sides q and b at each step
   *  \param time current time
   *  \return false if there is nothing to solve
   */
  bool preCompute(double time) override;

  /** write the numerics PlasticityProblem passed to plasticity_2d_driver (numerics .dat
   *  format, plasticity_printInFile) into <prefix>_<n>.dat, n being the number of the call
   *  of compute (1 for the first time step)
   *  \param prefix prefix of the file names (empty: no dump)
   *  \param steps the calls whose problem is written (all if empty)
   */
  void setProblemDump(const std::string& prefix, const std::vector<unsigned int>& steps = {}) {
    _dump_prefix = prefix;
    _dump_steps = steps;
  }

  /** compute the one-step problem (assembly by GlobalFrictionContact::preCompute,
   *  condensation, numerics solve, postCompute)
   *  \param time current time
   *  \return the numerics solver info
   */
  int compute(double time) override;
};

}  // namespace siconos::mechanics::fem::nonsmooth_formulations

#endif
