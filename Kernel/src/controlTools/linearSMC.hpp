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

/*! \file linearSMC.hpp
  \brief General interface to define an actuator
*/

#ifndef linearSMC_H
#define linearSMC_H

#include "SiconosKernel.hpp"
#include <boost/circular_buffer.hpp>

class linearSMC : public Actuator
{
private:
  /** default constructor */
  linearSMC();

  /** dimension of the state space */
  unsigned int _nDim;

  /** index for saving data */
  unsigned int _indx;

  /** control variable */
  SP::SimpleVector _u;

  /** the vector defining the surface (\f$ s = Cx\f$) */
  SP::SiconosMatrix _Csurface;

  /** the sensor that feed the controller */
  SP::controlSensor _sensor;

  /** the dynamical system we are controlling */
  SP::FirstOrderLinearDS _DS;

  /** boolean to determined if the controller has been correctly initialized */
  bool _initDone;

  /** current \f$\Deltat\f$ (or timeStep) */
  double _curDeltaT;

  /** t0 for the controller */
  double _t0;

  /** T final for the controller */
  double _T;

  /** XXX Some matrix */
  SP::SimpleMatrix _B;

  /** the Model for the controller */
  SP::Model _SMC;

  /** the dynamical system for the controller */
  SP::FirstOrderLinearDS _DS_SMC;

  /** the Relation for the Controller */
  SP::FirstOrderLinearTIR _relationSMC;

  /** the NonSmoothLaw for the controller */
  SP::NonSmoothLaw _nsLawSMC;

  /** */
  SP::Interaction _interactionSMC;

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

  /** easy access to lambda */
  SP::SiconosVector _lambda;

  /** easy access to lambda */
  SP::SiconosVector _xController;

  SP::SiconosVector sampledControl;

  SP::EventsManager _eventsManager;
public:

  /** Constructor with a TimeDiscretisation and a Model.
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param a SP::Model
   */
  linearSMC(int, SP::TimeDiscretisation, SP::Model);

  /** Constructor with a TimeDiscretisation, a Model and a set of Sensor.
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param a SP::Model
   * \param a set of Sensor linked to this Actuator.
   */
  linearSMC(int, SP::TimeDiscretisation, SP::Model, const Sensors&);

  /** destructor
   */
  virtual ~linearSMC();

  /** initialize actuator data.
   */
  virtual void initialize();

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * \f$u_k = u_{k-1} + c_1e_k + c_2e_{k-1} + c_3e_{k-2}\f$, where
   * \f{align*}c_1 &= K_P - \frac{K_D}/{\Delta t} + K_I \Delta t \\
   * c_2 &= -1 - \frac{2K_D}/{\Delta t} \\
   * c_3 &= \frac{K_D}/{\Delta t} \f{align*}
   */
  void actuate();

  /** Set the value of _Csurface to newValue
   * * \param a SimpleVector newValue
   */
  void setCsurface(const SiconosMatrix&);

  /** Set _Csurface to pointer newPtr
   * \param a SP::SimpleVector
   */
  void setCsurfacePtr(SP::SiconosMatrix);
};
DEFINE_SPTR(linearSMC)
#endif
