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

/*! \file Actuator.hpp
  \brief General interface to define an actuator
*/

#ifndef Actuator_H
#define Actuator_H

#include <string>
#include "RuntimeException.hpp"
#include "SiconosPointers.hpp"

#include "SiconosFwd.hpp"

#include "SiconosAlgebraTypeDef.hpp"

#include "ControlTypeDef.hpp"
#include "SiconosControlFwd.hpp"

/** Actuators Base Class

   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) February 09, 2007

   Abstract class, interface to user-defined actuators.

   An Actuator is dedicated to act on parameters of the Model
   (especially z param. in DynamicalSystem) according to some specific
   values recorded thanks to sensors. It gives an interface for User
   who can implement its own Actuator.  clearly define which data he
   needs to save.

   An Actuator handles a TimeDiscretisation, which defines the set of
   all instants where the Actuator must operate \n (i.e. each times
   where actuate() function will be called). An Event, inserted into
   the EventsManager of the Simulation, is linked to this
   TimeDiscretisation.

   Moreover, an Actuator is identified thanks to an id and a type (a
   number associated to the derived class type indeed).

   \section BActuator Construction

   To build an Actuator it is necessary to use the factory. Inputs are
   a number which identify the derived class type and a
   TimeDiscretisation:

   \code
   // Get the registry
   ActuatorFactory::Registry& regActuator(ActuatorFactory::Registry::get()) ;
   // Build an Actuator of type "myType" with t as a TimeDiscretisation.
   regActuator.instantiate(myType, t);
   \endcode

   The best way is to use the controlManager:
   \code
   // cm a ControlManager
   cm->addActuator(myType,t);
   // or if cm has already been initialized:
   cm->addAndRecordActuator(myType,t)
   \endcode
*/
class Actuator
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Actuator);

  /** type of the Actuator */
  unsigned int _type;

  /** id of the Actuator */
  std::string _id;

  /** Control variable */
  SP::SiconosVector  _u;

  /** B Matrix */
  SP::SimpleMatrix _B;

  /** ControlSensor feeding the Controller */
  SP::ControlSensor _sensor;

  /** default constructor
   */
  Actuator();

  /** copy constructor
   * Private => no copy nor pass-by value allowed.
   * \param a an Actuator
   */
  Actuator(const Actuator& a);

public:

  /** Constructor with a TimeDiscretisation.
   * \param type the type of the Actuator, which corresponds to the class type
   * \param sensor the ControlSensor feeding the Actuator
   */
  Actuator(unsigned int type, SP::ControlSensor sensor);

  /** destructor
   */
  virtual ~Actuator();

  /** set id of the Actuator
   *  \param newId the new id.
   */
  inline void setId(const std::string& newId)
  {
    _id = newId;
  };

  /** get id of the Actuator
   *  \return a std::string
   */
  inline const std::string getId() const
  {
    return _id;
  };

  /** get the type of the Actuator (ie class name)
   *  \return an integer
   */
  inline int getType() const
  {
    return _type;
  };

  /** Get the control value
   * \return current control value u
   */
  inline const SiconosVector& u() const { return *_u; };

  /** Set the B matrix
   * \param B the new B matrix
  */
  void setB(const SimpleMatrix& B);

  /** Set the B matrix
   * \param B the new B matrix
   */
  inline void setBPtr(SP::SimpleMatrix B)
  {
    _B = B;
  };

  /** add a Sensor in the actuator.
   *  \param newSensor a Sensor that will be connected to the Actuator
   */
  void addSensorPtr(SP::ControlSensor newSensor);

  /** This is derived in child classes if they need to copy the TimeDiscretisation
   * associated with this Actuator
  *  \param td the TimeDiscretisation for this Actuator
  */
  virtual void setTimeDiscretisation(const TimeDiscretisation& td) {};

  /** initialize actuator data.
   * \param m a SP::Model
   */
  virtual void initialize(const Model& m);

  /** capture data when the ActuatorEvent is processed
   */
  virtual void actuate() = 0;

  /** display the data of the Actuator on the standard output
   */
  virtual void display() const;

};
#endif
