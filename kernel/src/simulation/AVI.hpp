/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*! \file AVI.hpp
  Affine Variational Inequalities formulation
*/

#ifndef AVI_H
#define AVI_H

#include "AVI_cst.h"  // contains only enum. Ok.
#include "LinearOSNS.hpp"

struct AffineVariationalInequalities;

namespace siconos::nonsmooth_formulations {
/**
   Formalization and Resolution of an Affine Variational Inequality (AVI)

   This class is devoted to the formalization and the resolution of
   Affine variational Inequalities (AVI): given a polytopic set \f$ P \f$, \f$ M\in R^{p\times
   p} \f$ and \f$ q\in R^p \f$,

   \f[
   \text{find }z \in P\text{ such that}\quad \langle Mz+q, x - z\rangle \geq 0 \qquad \forall x
   \in P \f]

   \todo : add "recover" function to start from old values of z and w.

   For details regarding the available options, see Nonsmooth problems formulations and
   available solvers in users' guide.
*/

class AVI : public LinearOSNS {
 protected:
  ACCEPT_SERIALIZATION(AVI);

  /** contains the numerics problem for the AVI system */
  std::shared_ptr<AffineVariationalInequalities> _numerics_problem{nullptr};

 public:
  /** constructor from numerics solver id
   *
   *  \param numericsSolverId id of numerics solver, default =  SICONOS_AVI_CAOFERRIS
   */
  AVI(int numericsSolverId = SICONOS_AVI_CAOFERRIS);

  /** constructor from a pre-defined solver options set
   *
   *  \param options the options set
   *
   */
  AVI(std::shared_ptr<SolverOptions> options);

  /** destructor */
  virtual ~AVI() noexcept;

  void initialize(std::shared_ptr<siconos::simulation::Simulation> sim);

  /** Compute the unknown z and update the Interaction (y and lambda)
   *
   *  \param time current time
   *  \return information about the solver convergence.
   */
  int compute(double time);

  /** Check the compatibility fol the nslaw with the targeted OSNSP */
  bool checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw& nslaw);

  /** print the data to the screen
   */
  void display() const;
};
}  // namespace siconos::nonsmooth_formulations

#endif  // AVI_H
