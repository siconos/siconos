/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*! \file LinearSMCOT2.hpp
  \brief General interface to define a sliding mode controller with
  disturbance compensation. Reference: Su, W.C.; Drakunov, S.; Özgüner, Ü.
  An O(T2) boundary layer in sliding mode for sampled-data systems
  */

#ifndef LinearSMCOT2_H
#define LinearSMCOT2_H

#include "CommonSMC.hpp"
#include "OneStepIntegratorTypes.hpp"

class LinearSMCOT2 : public CommonSMC
{
private:
  /** default constructor */
  LinearSMCOT2() {};
  /** Current value of the state (\f$ x_k\f$)*/
  SP::SiconosVector _X;
  /** Predicted current value of the state (\f$ \hat{x}_k = \Phi x_{k-1} + \Gamma u_{k-1}\f$)*/
  SP::SiconosVector _Xhat;
  /** Next value of the state only with the influence of the dynamic \f$ \XPhi = \Phi x_k\f$*/
  SP::SiconosVector _XPhi;
  /** Model for the computation of _XPhi*/
  SP::NonSmoothDynamicalSystem _nsdsPhi;
  /** DynamicalSystem for the computation of _XPhi*/
  SP::FirstOrderLinearDS _DSPhi;
  /** TimeDiscretisation for the computation of _XPhi*/
  SP::TimeDiscretisation _tdPhi;
  /** OneSteoIntegrator for the computation of _XPhi*/
  SP::LsodarOSI _PhiOSI;
  /** Simulation for the computation of _XPhi*/
  SP::EventDriven _simulPhi;
  /** Model for the computation of Xhat*/
  SP::NonSmoothDynamicalSystem _nsdsPred;
  /** TimeDiscretisation for the computation of Xhat*/
  SP::TimeDiscretisation _tdPred;
  /** OneSteoIntegrator for the computation of Xhat*/
  SP::LsodarOSI _PredOSI;
  /** Simulation for the computation of Xhat*/
  SP::EventDriven _simulPred;
  /** DynamicalSystem for the computation of _Xhat*/
  SP::FirstOrderLinearDS _DSPred;
  /** Coefficient*/
  double _coeff;

  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LinearSMCOT2);

public:

  /** Constructor
   * \param sensor the ControlSensor feeding the Actuator
   */
  LinearSMCOT2(SP::ControlSensor sensor);

  /** destructor
  */
  virtual ~LinearSMCOT2();

  /** initialize actuator data.
   * \param m the Model
  */
  void initialize(const NonSmoothDynamicalSystem& nsds, const Simulation & s);

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * TODO
   */
  void actuate();

    /** This is derived in child classes if they need to copy the TimeDiscretisation
   * associated with this Sensor
  *  \param td the TimeDiscretisation for this Sensor
  */
  virtual void setTimeDiscretisation(const TimeDiscretisation& td);

};
#endif
