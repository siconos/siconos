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

#include <limits>

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

  enum MoreauJeanOSI_ds_workVector_id {
    RESIDU_FREE,
    VFREE,
    RESIDU_SIGMAFREE,
    SIGMAFREE,
    Q_SIGMAFREE,
    BUFFER,
    QTMP,
    WORK_LENGTH
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

  /** Apply the rule to one Interaction to know if it should be included in the
   *  IndexSet of level i
   *
   *  @param inter the Interaction to test
   *  @param i level of the IndexSet
   *  @return Boolean
   */
  virtual bool addInteractionInIndexSet(std::shared_ptr<siconos::modeling::Interaction> inter,
                                        unsigned int i) override;
};

}  // namespace siconos::mechanics::fem::integrators
#endif  // MoreauJeanOSI_H
