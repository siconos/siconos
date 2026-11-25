// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! @file  FEM-MoreauJeanGOSI.hpp
 * @brief specific implementation of the integrator for SolidLinearTIDS */

#ifndef FEMMoreauJeanGOSI_H
#define FEMMoreauJeanGOSI_H

#include "FEM-MoreauJeanOSI.hpp"

namespace siconos::mechanics::fem::integrators {

/** A specific implementation of MoreauJeanGOSI for SolidLinearTIDS */

class MoreauJeanGOSI : public MoreauJeanOSI {
 protected:
  ACCEPT_SERIALIZATION(MoreauJeanGOSI);

 public:
  using MoreauJeanOSI::MoreauJeanOSI;

  /** destructor */
  virtual ~MoreauJeanGOSI() noexcept = default;

  siconos::algebra::SiconosVector& get_v_iter(
      std::vector<std::shared_ptr<algebra::SiconosVector>> ds_works) {
    return *ds_works[tools::enum_to_index(wk_ds::v_iter)];
  }
  siconos::algebra::SiconosVector& get_sigma_iter(
      std::vector<std::shared_ptr<algebra::SiconosVector>> ds_works) {
    return *ds_works[tools::enum_to_index(wk_ds::sigma_iter)];
  }

  /** \return the maximum of all norms for the "MoreauJeanGOSI-discretized" residus of DS
   */
  virtual double computeResidu() override;

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input (_r or _p)
   */
  void computeFreeState() override {};

  /** integrate the system, between tinit and tend (->iout=true), with possible
   *  stop at tout (->iout=false)
   *
   *  \param tinit the initial time
   *  \param tend the end time
   *  \param tout the real end time
   *  \param notUsed useless flag (for MoreauJeanGOSI, used in LsodarOSI)
   */
  void integrate(double& tinit, double& tend, double& tout, int& notUsed) override;

  /** compute the current iteration
   */
  void computeIteration() override;

  /** Compute the nonsmooth law contribution to the output
   *
   *  \param inter the interaction (for y_k)
   *  \param osnsp the non-smooth integrator
   */
  void NonSmoothLawContributionToOutput(
      std::shared_ptr<siconos::modeling::Interaction> inter,
      siconos::nonsmooth_formulations::OneStepNSProblem& osnsp);
};
}  // namespace siconos::mechanics::fem::integrators
#endif  // MoreauJeanGOSI_H
