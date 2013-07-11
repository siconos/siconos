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

/*! \file CommonSMC.hpp
  \brief General interface to define a sliding mode actuator
  */

#ifndef CommonSMC_H
#define CommonSMC_H

#include "Actuator.hpp"
#include "SiconosAlgebraTypeDef.hpp"

#ifndef FirstOrderLinearDS_H
DEFINE_SPTR(FirstOrderLinearDS)
#endif
#ifndef TimeStepping_H
DEFINE_SPTR(TimeStepping)
#endif
#ifndef Relay_H
DEFINE_SPTR(Relay)
#endif
#ifndef EventsManager_H
DEFINE_SPTR(EventsManager)
#endif

class CommonSMC : public Actuator
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(CommonSMC);

protected:
  /** default constructor */
  CommonSMC() {};

  /** index for saving data */
  unsigned int _indx;

  /** the vector defining the surface (\f$ s = Cx \f$) */
  SP::SiconosMatrix _Csurface;

  /** matrix describing the influence of \f$lambda\f$ on s */
  SP::SiconosMatrix _D;

  /** scalar multiplying Sign; \f$ u^s = - \alpha Sign */
  double _alpha;

  /** the Relation for the Controller */
  SP::Relation _relationSMC;

  /** the NonSmoothLaw for the controller */
  SP::NonSmoothLaw _sign;

  /** Interaction for the control */
  SP::Interaction _interactionSMC;

  /** easy access to lambda */
  SP::SiconosVector _lambda;

  /** Relay solver type */
  int _numericsSolverId;

  /** Numerical precision expected for the Relay solver */
  double _precision;

  /** the Model for the controller */
  SP::Model _SMC;
  /** the DynamicalSystem for the controller */
  SP::FirstOrderLinearDS _DS_SMC;
  /** the TimeDiscretisation for the controller */
  SP::TimeDiscretisation _td;
  /** Simulation for the controller */
  SP::TimeStepping _simulationSMC;
  /** Integrator for the controller */
  SP::OneStepIntegrator _integratorSMC;
  /** Theta for the controller */
  double _thetaSMC;
  /** OneStepNsProblem for the controller */
  SP::Relay _OSNSPB_SMC;

  /** SP::EventsManager of the SMC Simulation */
  SP::EventsManager _eventsManager;
  /** SP::NonSmoothLaw for computing the control law */
  SP::NonSmoothLaw _nsLawSMC;

  /** inverse of CB */
  SP::SimpleMatrix _invCB;

  /** Store \f$u^{eq}\f$ */
  SP::SiconosVector _ueq;

  /** Store \f$u^s\f$ */
  SP::SiconosVector _us;

  /** Do not use the state-continuous equivaluent control \f$u^{eq}\f$ */
  bool _noUeq;

  /** Compute the equivalent part of the control \f$u^{eq}\f$.
   * The method used here is to discretize the continuous-time
   * formula using a theta method
   */
  void computeUeq();

public:

  /** Constructor with a TimeDiscretisation and a Model.
   * \param type the type of the SMC Actuator
   * \param sensor the ControlSensor feeding the Actuator
   */
  CommonSMC(unsigned int type, SP::ControlSensor sensor): Actuator(type, sensor),
    _alpha(1.0), _numericsSolverId(SICONOS_RELAY_LEMKE), _precision(1e-8),
    _thetaSMC(0.5), _noUeq(false) {}

  /** Constructor with a TimeDiscretisation, a Model and two matrices
   * \param type the type of the SMC Actuator
   * \param sensor the ControlSensor feeding the Actuator
   * \param B the B matrix
   * \param D the saturation matrix
   */
  CommonSMC(unsigned int type, SP::ControlSensor sensor, SP::SiconosMatrix B, SP::SiconosMatrix D):
    Actuator(type, sensor), _D(D), _alpha(1.0), _numericsSolverId(SICONOS_RELAY_LEMKE),
    _precision(1e-8), _thetaSMC(0.5), _noUeq(false)
  {
    _B = B;
  }


  /** Compute the new control law at each event
  */
  virtual void actuate() = 0;

  /** Initialization
   * \param m a SP::Model
   */
  virtual void initialize(const Model& m);

  /** Set the value of _Csurface to newValue
   * \param newValue the new value for _Csurface
   */
  void setCsurface(const SiconosMatrix& newValue);

  /** Set _Csurface to pointer newPtr
   * \param newPtr a SP::SiconosMatrix containing the new value for _Csurface
   */
  void setCsurfacePtr(SP::SiconosMatrix newPtr);

  /** Set the value of _D to newValue
   * \param newValue the new value for _D
   */
  void setSaturationMatrix(const SiconosMatrix& newValue);

  /** Set _D to pointer newPtr
   * \param newPtr a SP::SiconosMatrix containing the new value for _D
   */
  void setSaturationMatrixPtr(SP::SiconosMatrix newPtr);

  /** Set _alpha
   * \param alpha the new value for _alpha
   */
  inline void setAlpha(double alpha) { _alpha = alpha; };

  /** get _lambda
   * \return a pointer to _lambda
   */
  inline SP::SiconosVector lambda() const
  {
    return _lambda;
  };

  /** Set the solver
   * \param newNumericsSolverId the solver for the relay
   */
  inline void setSolver(const int newNumericsSolverId)
  {
    _numericsSolverId = newNumericsSolverId;
  };

  /** Set the precision
   * \param newPrecision a double
   */
  inline void setPrecision(const double newPrecision)
  {
    _precision = newPrecision;
  };

  /** Get the Relay problem associated with the controller. This is useful to
   * gain access to the data given to the Relay solver in Numerics
   * \return a reference to the Relay problem
   */
  inline const Relay& relay()
  {
    return * _OSNSPB_SMC;
  };

  /** get \f$u^{eq}\f$
   * \return a reference to _ueq
   */
  inline SiconosVector& ueq()
  {
    return *_ueq;
  };

  /** get \f$u^{s}\f$
   * \return a reference to _us
   */

  inline SiconosVector& us()
  {
    return *_us;
  };

  /** Set _theta, used in some discretization method for \f$u^{eq}\f$
   * \param newTheta the new value for _thetaSMC
   */

  inline void setTheta(const double& newTheta)
  {
    _thetaSMC = newTheta;
  };

  /** Disable (or enable) the use of the state-continuous control \f$u^{eq}\f$
   * \param b disable the use of Ueq if true
   */
  inline void noUeq(const bool b)
  {
    _noUeq = b;
  };

  /** This is derived in child classes if they need to copy the TimeDiscretisation
   * associated with this Sensor
  *  \param td the TimeDiscretisation for this Sensor
  */
  virtual void setTimeDiscretisation(const TimeDiscretisation& td);

};
DEFINE_SPTR(CommonSMC)
#endif
