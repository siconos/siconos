/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

/*! \file LinearSMC.hpp
  \brief General interface to define an actuator
*/

#ifndef LinearSMC_H
#define LinearSMC_H

#include "SiconosKernel.hpp"
#include <boost/circular_buffer.hpp>

class LinearSMC : public CommonSMC
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(LinearSMC);

  /** default constructor */
  LinearSMC() {};

  /** B Matrix in the Interaction */
  SP::SiconosMatrix _B;
  /** D Matrix in the Interaction */
  SP::SiconosMatrix _D;
  /** the Model for the controller */
  SP::Model _SMC;
  /** the DynamicalSystem for the controller */
  SP::FirstOrderLinearDS _DS_SMC;
  /** the TimeDiscretisation for the controller */
  SP::TimeDiscretisation _tD_SMC;
  /** Simulation for the controller */
  SP::TimeStepping _simulationSMC;
  /** Integrator for the controller */
  SP::Moreau _integratorSMC;
  /** Theta for the controller */
  double _thetaSMC;
  /** LCP for the controller */
  SP::LCP _LCP_SMC;
  /** OneStepNsProblem for the controller */
  SP::Relay _OSNSPB_SMC;
  /** SP::SiconosVector containing the control */
  SP::SiconosVector _sampledControl;
  /** SP::EventsManager of the SMC Simulation */
  SP::EventsManager _eventsManager;
  /** SP::NonSmoothLaw for computing the control law */
  SP::NonSmoothLaw _nsLawSMC;

public:

  /** Constructor with a TimeDiscretisation and a DynamicalSystem.
   * \param t a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param ds the SP::DynamicalSystem we are controlling
   */
  LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds);

  /** Constructor with a TimeDiscretisation and a DynamicalSystem.
   * \param t a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param ds the SP::DynamicalSystem we are controlling
   * \param B the B matrix in the FirstOrderLinearR
   * \param D the D matrix in the FirstOrderLinearR
   */
  LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, SP::SiconosMatrix B, SP::SiconosMatrix D);

  /** Constructor with a TimeDiscretisation, a DynamicalSystem and a set of Sensor.
   * \param t a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param ds the SP::DynamicalSystem we are controlling
   * \param sensorList a set of Sensor linked to this Actuator.
   */
  LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList);

  /** Initialize
   * \param m a SP::Model
   */
  void initialize(SP::Model m);

  /** destructor
   */
  virtual ~LinearSMC();

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * TODO
   */
  void actuate();

  /** Set the B matrix
   * \param B the new B matrix
  */
  inline void setB(SiconosMatrix & B)
  {
    _B.reset(new SimpleMatrix(B));
  };
  /** Set the B matrix
   * \param B the new B matrix
  */
  inline void setBPtr(SP::SiconosMatrix B)
  {
    _B = B;
  };

  /** Set the D matrix
   * \param D the new D matrix
  */
  inline void setD(SiconosMatrix & D)
  {
    _D.reset(new SimpleMatrix(D));
  };
  /** Set the D matrix
   * \param D the new D matrix
  */
  inline void setDPtr(SP::SiconosMatrix D)
  {
    _D = D;
  };


};
DEFINE_SPTR(LinearSMC)
#endif
