/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

#ifndef FirstOrderLinearDS_H
DEFINE_SPTR(FirstOrderLinearDS)
#endif
#ifndef EventDriven_H
DEFINE_SPTR(EventDriven)
#endif
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
  SP::Model _modelPhi;
  /** DynamicalSystem for the computation of _XPhi*/
  SP::FirstOrderLinearDS _DSPhi;
  /** TimeDiscretisation for the computation of _XPhi*/
  SP::TimeDiscretisation _tdPhi;
  /** OneSteoIntegrator for the computation of _XPhi*/
  SP::Lsodar _PhiOSI;
  /** Simulation for the computation of _XPhi*/
  SP::EventDriven _simulPhi;
  /** Model for the computation of Xhat*/
  SP::Model _modelPred;
  /** TimeDiscretisation for the computation of Xhat*/
  SP::TimeDiscretisation _tdPred;
  /** OneSteoIntegrator for the computation of Xhat*/
  SP::Lsodar _PredOSI;
  /** Simulation for the computation of Xhat*/
  SP::EventDriven _simulPred;
  /** DynamicalSystem for the computation of _Xhat*/
  SP::FirstOrderLinearDS _DSPred;
  /** Coefficient*/
  double _coeff;

protected:
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
  void initialize(const Model& m);

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * TODO
   */
  void actuate();

    /** This is derived in child classes if they need to copy the TimeDiscretisation
   * associated with this Sensor
  *  \param td the TimeDiscretisation for this Sensor
  */
  virtual void setTimeDiscretisation(const TimeDiscretisation& td)
  { 
    _tdPhi.reset(new TimeDiscretisation(td));
    _tdPred.reset(new TimeDiscretisation(td));
  };

};
DEFINE_SPTR(LinearSMCOT2)
#endif
