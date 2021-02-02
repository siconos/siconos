/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

/*! \file Actuator.hpp
  \brief General interface to define an actuator
*/

#ifndef Actuator_H
#define Actuator_H

#include <string>
#include "SiconosPointers.hpp"

#include "SiconosFwd.hpp"

#include "SiconosAlgebraTypeDef.hpp"

#include "ControlTypeDef.hpp"
#include "SiconosControlFwd.hpp"

/** Actuators Base Class

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

private:

  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Actuator);

protected:

  /** type of the Actuator */
  unsigned int _type;

  /** id of the Actuator */
  std::string _id;

  /** Control variable */
  SP::SiconosVector  _u;

  /** B Matrix */
  SP::SimpleMatrix _B;

  /** name of the plugin for g (nonlinear affine in control system)*/
  std::string _plugingName;

  /** name of the plugin to compute \f$\nabla_x g\f$ for the nonlinear case*/
  std::string _pluginJacgxName;

  /** ControlSensor feeding the Controller */
  SP::ControlSensor _sensor;

  /** default constructor
   */
  Actuator();

public:

  /** General Constructor
   * \param type the type of the Actuator, which corresponds to the class type
   * \param sensor the ControlSensor feeding the Actuator
   */
  Actuator(unsigned int type, SP::ControlSensor sensor);

  /** General Constructor with dynamics affine in control
   * \param type the type of the Actuator, which corresponds to the class type
   * \param sensor the ControlSensor feeding the Actuator
   */
  Actuator(unsigned int type, SP::ControlSensor sensor, SP::SimpleMatrix B);

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

  /** Set the control size
   * \param size dimension of the control input u
   */
  void setSizeu(unsigned size);

  /** Set the B matrix
   * \param B the new B matrix
   */
  inline void setB(SP::SimpleMatrix B)
  {
    _B = B;
  };

  /** Set the name of the plugin for computing g
   * \param g the name of the plugin to compute g
   */
  inline void setg(const std::string& g)
  {
    _plugingName = g;
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
   * \param nsds the NonSmoothDynamicalSystem
   * \param s the simulation
   */
  virtual void initialize(const NonSmoothDynamicalSystem& nsds, const Simulation & s);

  /** capture data when the ActuatorEvent is processed
   */
  virtual void actuate() = 0;

  /** display the data of the Actuator on the standard output
   */
  virtual void display() const;

  /** get the NSDS used in the Controller, if there is one
   * \return "NULL" shared_ptr if there is no internal simulation, otherwise
   * it return the Model hoding the simulation
   */
  virtual SP::NonSmoothDynamicalSystem getInternalNSDS() const;
};
#endif
