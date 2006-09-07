/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#ifndef EVENT_H
#define EVENT_H

/** \class Event
 *  \brief virtual class that represents generic time events.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date (Creation) February 21, 2006
 *
 *  This base class simply records the time at which the event will take place. A pure virtual function named process
 *  will be invoked to execute the event.
 *  The time is represented with an unsigned long int.
 *  TODO: define with more details what process is supposed to do in any case.
 */

#include "SiconosConst.h"
#include "RuntimeException.h"
#include<string>
#include<iostream>

class Simulation;

const unsigned long int DEFAULT_EVENT_TIME = 0;
const std::string DEFAULT_EVENT_TYPE = "undefined";

class Event
{
protected:

  /** Date of the present event,
   *  represented with a long int */
  const unsigned long int timeOfEvent;

  /** Id or type of the Event */
  const std::string type;

  /** Default constructor */
  Event();

public:

  /** \fn Event(const Event&)
   *  \brief copy constructor
   *  \param the event to be copied
   */
  Event(const Event&);

  /** \fn Event(const unsigned long int&, const string &)
   *  \brief constructor with time value and type as input
   *  \param an unsigned int
   *  \param a string
   */
  Event(const unsigned long int&, const std::string & = DEFAULT_EVENT_TYPE);

  /** \fn ~Event()
   *  \brief destructor
   */
  virtual ~Event();

  // GETTERS/SETTERS

  /** \fn const unsigned long int getTimeOfEvent() const
   *  \brief get the time of the present event (unsigned int format)
   *  \return an unsigned long int
   */
  inline const unsigned long int getTimeOfEvent() const
  {
    return timeOfEvent ;
  };

  /** \fn void setIntTimeOfEvent(const unsigned long int&)
   *  \brief set the timeOfEvent as a long int
   *  \param an unsigned long int
   */
  //  inline void setIntTimeOfEvent(const unsigned long int newVal) {timeOfEvent = newVal;};

  /** \fn const std::string getType() const
   *  \brief get a type of the present event
   *  \return an std::string
   */
  inline const std::string getType() const
  {
    return type ;
  };

  /** \fn void setType(const std::string&)
   *  \brief set the type of this event
   *  \param an std::string
   */
  //inline void setType(const std::string newVal) {type = newVal;};

  /** \fn const double getTimeOfEvent() const
   *  \brief get the time (double format) of present event
   *  \return a double
   */
  //const double getTimeOfEvent() const;

  /** \fn void setTimeOfEvent(const double&)
   *  \brief set dateOfEvent using a double input for time
   *  \param a double
   */
  //void setTimeOfEvent(const double& newVal);

  /** \fn void saveEventToXML()
   *  \brief copy the data of the Event into the XML tree
   */
  // void saveEventToXML();

  /** \fn void display()
   *  \brief display Event data
   */
  void display() const ;

  /** \fn virtual void process(Simulation*) = 0
   *  \brief virtual function which actions depends on event type
   * \param Simulation*, the simulation that owns this Event (through the EventsManager)
   */
  virtual void process(Simulation*) = 0;
};

#endif // Event_H
