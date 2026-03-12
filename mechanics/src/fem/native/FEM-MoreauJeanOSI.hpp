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

/*! \file  FEM-MoreauJeanOSI.hpp */

#ifndef FEMMoreauJeanOSI_H
#define FEMMoreauJeanOSI_H

#include "MoreauJeanOSI.hpp"

namespace siconos::mechanics::fem::integrators {

/**
   One Step time Integrator, Moreau-Jean algorithm dedicated to nsds based
   on DS and relations from the mechanics::fem component.

   For details regarding MoreauJeanOSI, see siconos::integrators::MoreauJeanOSI

*/
class MoreauJeanOSI : public siconos::integrators::MoreauJeanOSI {
 protected:
  ACCEPT_SERIALIZATION(MoreauJeanOSI);

 public:
  using siconos::integrators::MoreauJeanOSI::MoreauJeanOSI;

  /** This enum is used to get access to work vectors relared to DS
   *  It corresponds to:
   *  - a container saved in the graph of ds
   *  - a container associated to a specific ds
   */
  enum class wk_ds : std::size_t {
    residu_free,
    vfree,
    v_iter,
    residu_sigma_free,
    sigma_free,
    sigma_iter,
    q_sigma_free,
    buffer,
    size
  };

  /** destructor */
  virtual ~MoreauJeanOSI() noexcept = default;

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   *
   *  @param t time of initialization
   *  @param ds the dynamical system
   */
  virtual void initializeWorkVectorsForDS(
      double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) override;

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

  /** initialize iteration matrix W MoreauJeanOSI matrix at time t
   *
   *  @param time
   *  @param ds a pointer to DynamicalSystem
   */
  virtual void initializeIterationMatrix(
      double time, std::shared_ptr<siconos::modeling::SecondOrderDS> ds) const override;

  /** @return the inverse of the iteration matrix
   *  @param ds the concerned dynamical system
   *  @param LUW the LU-factorisation of W
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> iterationMatrixInverse(
      std::shared_ptr<siconos::modeling::SecondOrderDS> ds,
      const siconos::algebra::SiconosDenseLUMatrix& LUW);

  /** compute the initial state of the Newton loop */
  virtual void computeInitialNewtonState() override;

  /**
      return the maximum of all norms for the "MoreauJeanOSI-discretized"
      residus of DS

      @return a double
    */
  double computeResidu() override;

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input (_r or _p)
   */
  virtual void computeFreeState() override;

  /** integrate the system, between tinit and tend (->iout=true), with possible
   *  stop at tout (->iout=false)
   *
   *  @param tinit the initial time
   *  @param tend the end time
   *  @param tout the real end time
   *  @param notUsed useless flag (for MoreauJeanOSI, used in LsodarOSI)
   */
  void integrate(double& tinit, double& tend, double& tout, int& notUsed) override;

  /** update the state of the dynamical systems
   *
   *  @param level the level of interest for the dynamics: not used at the time
   */
  virtual void updateState(const unsigned int level) override;

  /** compute the current iteration
   *
   */
  void computeIteration() override;

  /** Apply the rule to one Interaction to know if it should be included in the
   *  IndexSet of level i
   *
   *  @param inter the Interaction to test
   *  @param i level of the IndexSet
   *  @return Boolean
   */
  virtual bool addInteractionInIndexSet(
      std::shared_ptr<siconos::modeling::Interaction> inter,
      siconos::graphs::InteractionsGraph::size_type i) override;
};

}  // namespace siconos::mechanics::fem::integrators
#endif  // MoreauJeanOSI_H
