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

/*! @file  FEM-MoreauJeanGOSI.hpp
 * @brief specific implementation of the integrator for SolidLinearTIDS */

#ifndef FEMMoreauJeanGOSI_H
#define GEMMoreauJeanGOSI_H

#include "MoreauJeanGOSI.hpp"

namespace siconos::mechanics::fem::integrators {

/** A specific implementation of MoreauJeanGOSI for SolidLinearTIDS */

class MoreauJeanGOSI : public siconos::integrators::MoreauJeanGOSI {
 protected:
  ACCEPT_SERIALIZATION(MoreauJeanGOSI);

 public:
  // Warning: enum could be mixed up with those of MoreauJeanOSI
  enum MoreauJeanGOSI_ds_workVector_id {
    RESIDU_FREE,
    FREE,
    RESIDU_SIGMAFREE,
    SIGMAFREE,
    Q_SIGMAFREE,
    LOCAL_BUFFER,
    WORK_LENGTH
  };

  using siconos::integrators::MoreauJeanGOSI::MoreauJeanGOSI;

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

  /** \return the maximum of all norms for the "MoreauJeanGOSI-discretized" residus of DS
   */
  virtual double computeResidu() override;
};
}  // namespace siconos::mechanics::fem::integrators
#endif  // MoreauJeanGOSI_H
