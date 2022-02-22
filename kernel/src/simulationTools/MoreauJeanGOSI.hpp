/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "MoreauJeanOSI.hpp"
#include "OneStepNSProblem.hpp"

#include <limits>

/**  \class MoreauJeanGOSI
 *   \brief A global version of the MoreauJeanOSI integrator
 */

class MoreauJeanGOSI : public MoreauJeanOSI {
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(MoreauJeanGOSI);

public:
  // Warning: enum could be mixed up with those of MoreauJeanOSI
  enum MoreauJeanGOSI_ds_workVector_id {
    RESIDU_FREE,
    FREE,
    LOCAL_BUFFER,
    WORK_LENGTH
  };

  // enum MoreauJeanGOSI_interaction_workVector_id{OSNSP_RHS,
  // WORK_INTERACTION_LENGTH};

  /** constructor from theta value only
   *  \param theta value for all linked DS (default = 0.5).
   *  \param gamma value for all linked DS (default = NaN and gamma is not
   * used).
   */
  MoreauJeanGOSI(double theta = 0.5,
                 double gamma = std::numeric_limits<double>::quiet_NaN())
      : MoreauJeanOSI(theta, gamma){};

  /** destructor
   */
  virtual ~MoreauJeanGOSI(){};

  // --- OTHER FUNCTIONS ---

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   * \param t time of initialization
   * \param ds the dynamical system
   */
  void initializeWorkVectorsForDS(double t, SP::DynamicalSystem ds) override;

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   * \param inter the interaction
   * \param interProp the properties on the graph
   * \param DSG the dynamical systems graph
   */
  void initializeWorkVectorsForInteraction(Interaction &inter,
                                           InteractionProperties &interProp,
                                           DynamicalSystemsGraph &DSG) override;

  /** return the maximum of all norms for the "MoreauJeanGOSI-discretized"
     residus of DS \return a double
   */
  double computeResidu() override;

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input (_r or _p)
   */
  void computeFreeState() override;

  /** integrate the system, between tinit and tend (->iout=true), with possible
   * stop at tout (->iout=false) \param tinit the initial time \param tend the
   * end time \param tout the real end time \param notUsed useless flag (for
   * MoreauJeanGOSI, used in LsodarOSI)
   */
  void integrate(double &tinit, double &tend, double &tout,
                 int &notUsed) override;

  /** update the state of the dynamical systems
   *  \param level the level of interest for the dynamics: not used at the time
   */
  void updateState(const unsigned int level) override;

  /** Compute the nonsmooth law contribution to the output
   * \param inter the interaction (for y_k)
   * \param osnsp the non-smooth integrator
   */
  void NonSmoothLawContributionToOutput(SP::Interaction inter,
                                        OneStepNSProblem &osnsp);

  /** Displays the data of the MoreauJeanGOSI's integrator
   */
  void display() override;

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();
};

#endif // MoreauJeanGOSI_H
