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

/*! \file Observer.hpp
  \brief General interface to define an Observer
*/

#ifndef Observer_H
#define Observer_H

#include <string>
#include <set>
#include "RuntimeException.hpp"
#include "SiconosPointers.hpp"

#include "SiconosFwd.hpp"

#include "SiconosAlgebraTypeDef.hpp"
#include "ControlTypeDef.hpp"

#ifndef TimeStepping_H
DEFINE_SPTR(TimeStepping)
#endif

/** Observers Base Class

   \author SICONOS Development Team - copyright INRIA
   \version 3.6.0
   \date (Creation) May 21, 2013

   Abstract class, interface to user-defined observers.

   An Observer is dedicated to estimate the state of a DynamicalSystem given
   its dynamics, inputs and a initial estimate of the state.
*/

class Observer : public std11::enable_shared_from_this<Observer>
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Observer);

  /** type of the Observer */
  int _type;

  /** the DynamicalSystem used in the Observer */
  SP::DynamicalSystem _DS;

  /** The TimeDiscretisation */
  SP::TimeDiscretisation _td;

  /** the sensor that feed the observer */
  SP::ControlSensor _sensor;

  /** estimated state */
  SP::SiconosVector _xHat;

  /** The error \f$e=\hat{y}-y\f$ */
  SP::SiconosVector _e;

  /** The measurements from the sensor */
  SP::SiconosVector _y;

  /** id of the Observer */
  std::string _id;

  /** Model for integration */
  SP::Model _model;

  /** Simulation for integration */
  SP::TimeStepping _simulation;

  /** Integration for integration */
  SP::OneStepIntegrator _integrator;

  /** default constructor
   */
  Observer();

public:

  /** Constructor with a TimeDiscretisation.
   * \param type the type of the Observer, which corresponds to the class type
   * \param sensor the SP::Sensor to get the measurements
   * \param xHat0 the initial guess for the state
   * \param newId the id of the Observer
   */
  Observer(unsigned int type, SP::ControlSensor sensor, const SiconosVector& xHat0, const std::string& newId = "none");

  /** Constructor with a TimeDiscretisation.
   * \param type the type of the Observer, which corresponds to the class type.
   * \param sensor the SP::Sensor to get the measurements
   * \param xHat0 the initial guess for the state
   * \param ds the SP::DynamicalSystem used as a model for the real DynamicalSystem
   * \param newId the id of the Observer
   */
  Observer(unsigned int type, SP::ControlSensor sensor, const SiconosVector& xHat0, SP::DynamicalSystem ds, const std::string& newId = "none");

  /** destructor
   */
  virtual ~Observer();

  /** set id of the Observer
   *  \param newId the new id.
   */
  inline void setId(const std::string& newId)
  {
    _id = newId;
  };

  /** get id of the Observer
   *  \return a string
   */
  inline const std::string getId() const
  {
    return _id;
  };

  /** get the type of the Observer (ie class name)
   *  \return an integer
   */
  inline int getType() const
  {
    return _type;
  };

  /** This is derived in child classes if they need to copy the TimeDiscretisation
   * associated with this Sensor
  *  \param td the TimeDiscretisation for this Sensor
  */
  virtual void setTimeDiscretisation(const TimeDiscretisation& td);

  /** initialize observer data.
   * \param m a SP::Model
   */
  virtual void initialize(const Model& m);

  /** capture data when the ObserverEvent is processed
   */
  virtual void process() = 0;

  /** display the data of the Observer on the standard output
   */
  void display() const;

  /** get the error e
   * \return a pointer to e
   */
  inline SP::SiconosVector e()
  {
    return _e;
  }

  /** get the estimated state
   * \return a pointer to xHat
   */
  inline SP::SiconosVector xHat()
  {
    return _xHat;
  }

  /** Set the DynamicalSystem used in the Observer
   * \param ds the DynamicalSystem used in the Observer
   */
  inline void setDS(SP::DynamicalSystem ds)
  {
    _DS = ds;
  }

};

#endif
