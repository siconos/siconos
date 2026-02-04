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

/*! \file  MoreauJeanGOSI.hpp */

#ifndef MoreauJeanGOSI_H
#define MoreauJeanGOSI_H

#include <limits>
#include <memory>

#include "MoreauJeanOSI.hpp"
#include "SiconosVector.hpp"

namespace siconos::integrators {

/**
    A global version of the MoreauJeanOSI integrator
*/

class MoreauJeanGOSI : public MoreauJeanOSI {
 protected:
  ACCEPT_SERIALIZATION(MoreauJeanGOSI);

 public:
  /** constructor from theta value only
   *
   *  \param theta value for all linked DS (default = 0.5).
   *  \param gamma value for all linked DS (default = NaN and gamma is not
   *  used).
   */
  MoreauJeanGOSI(double theta = 0.5, double gamma = std::numeric_limits<double>::quiet_NaN())
      : MoreauJeanOSI(theta, gamma) {};

  /** destructor */
  virtual ~MoreauJeanGOSI() noexcept = default;

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   *
   *  \param t time of initialization
   *  \param ds the dynamical system
   */
  virtual void initializeWorkVectorsForDS(
      double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) override;

  virtual siconos::algebra::SiconosVector& get_v_iter(
      std::vector<std::shared_ptr<algebra::SiconosVector>> ds_works) {
    return *ds_works[tools::enum_to_index(wk_ds::v_iter)];
  };

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   *
   *  \param inter the interaction
   *  \param interProp the properties on the graph
   *  \param DSG the dynamical systems graph
   */
  virtual void initializeWorkVectorsForInteraction(
      siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
      siconos::graphs::DynamicalSystemsGraph& DSG) override;

  /** compute the initial state of the Newton loop.
   */
  void computeInitialNewtonState() override;

  /** \return the maximum of all norms for the "MoreauJeanGOSI-discretized" residus of DS
   */
  virtual double computeResidu() override;

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input (_r or _p)
   */
  void computeFreeState() override;

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

  /** Displays the data of the MoreauJeanGOSI's integrator
   */
  void display() const override;
};
}  // namespace siconos::integrators
#endif  // MoreauJeanGOSI_H
