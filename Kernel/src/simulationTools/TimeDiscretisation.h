/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
/*! \file TimeDiscretisation.h

  \brief A time-discretisation scheme for a given interval.
 */
#ifndef TIMEDISCRETISATION_H
#define TIMEDISCRETISATION_H

#include "TimeDiscretisationXML.h"
#include <vector>

class Model;
class TimeDiscretisationXML;

/** A time discretisation scheme

    \author SICONOS Development Team - copyright INRIA
    \version 3.0.0.
    \date (Creation) Apr 26, 2004

    A TimeDiscretisation object is used to discretized a given time interval. \n
    TimeDiscretisation are used:
    - in the simulation, as a user-input to discretized [t0,T]
    - in Sensor and Actuator, to define the set of time instants where the sensor or actuator
    must operate.

    A TimeDiscretisation is defined with a starting time (t0), a time step size (h, non necessarily constant),
    the number of the current time step (k). \n
    The time instant values are saved in a vector tk. Depending on the way of construction of the TimeDiscretisation,
    all or only current and next times are saved in tk. The difference is just a question of saving memory. \n

    Note that the TimeDiscretisation is linked to a Model. This is useful to synchronise it with initial and final
    times of the Model.

    \section tdMfunc Main functions:
    - setCurrentTimeStep(), to set current h. This value will be used for all future time steps, until next change.
    - increment(), shift to next time step (increment k, and shift t[k] and t[k+1])
    - getCurrentTime(), return t[k]

    \section tdBuild Construction

    - input = the complete vector tk. This defines t0, T, number of time steps and time step size
    (which is not necessarily constant). In this case, the whole vector is saved in the memory.
    - inputs = number of time steps and the Model. t0 and T are given by the input model, and the time step
    size h is computed with t0,T and nSteps. Only two values are saved: t[k] and t[k+1] = t[k] + h.
    h can be changed at any time.
    - inputs = h and the Model. t0 is given by the input model. Only two values are saved: t[k] and t[k+1] = t[k] + h.
    h can be changed at any time.
    - inputs = t0 and h.  Only two values are saved: t[k] and t[k+1] = t[k] + h.
    h can be changed at any time.

 */
class TimeDiscretisation
{
private:

  /** Default value for the time step (tk+1 - tk) */
  double h;

  /** Number of current time step (Simulation runs between steps k and k+1) */
  unsigned int k;

  /** vector of time values at each step (=> size = n1+n2+1 - Default size = 2 - Max size= nSteps+1) */
  TkVector tk;

  /** the XML object linked to the TimeDiscretisation to read XML data */
  TimeDiscretisationXML* timeDiscretisationXML;

  /* the model linked to this discretisation (used to set t0 and possibly final T) */
  Model* model;

  /** Flag to check if all data are up to date (ie if it is necessary to initialize or not) */
  bool isUpToDate;

  /** Indic. flag which sets the way the time-discretisation is built.*/
  int tdCase;

  /** index in tk which corresponds to the current time step (tk[pos] = t[k]) -
      Required since the size of tk depends on tdCase.
      tdCase = 1 => pos = k - tdCase = 2 => pos = 0 .
  */
  int pos;

  /** default constructor (private => no copy nor pass-by value)
   */
  TimeDiscretisation();

  /** Assignment Operator (private => forbidden) */
  TimeDiscretisation& operator =(const TimeDiscretisation&);

public:

  /** constructor with XML
   *  \param TimeDiscretisationXML* : the XML object corresponding
   *  \param Model* : the model that owns this discretisation
   */
  TimeDiscretisation(TimeDiscretisationXML*, Model *);

  // --- Straightforward constructors ---

  /** constructor with tk and model as given data
   *  \param a TkVector that describes the discretisation
   *  \param Model* : the model that owns this discretisation
   */
  TimeDiscretisation(const TkVector&, Model*);

  /** constructor with the number of time steps
   *  \param int (nb steps)
   *  \param Model* : the model that owns this discretisation
   */
  TimeDiscretisation(unsigned int, Model*);

  /** constructor with the size of the default time step
   *  \param double, the time step
   *  \param Model* : the model that owns this discretisation
   */
  TimeDiscretisation(double, Model*);

  /** constructor with the size of the default time step and t0
   *  \param double, initial time value
   *  \param double, the time step
   */
  TimeDiscretisation(double, double);

  // Destructor
  ~TimeDiscretisation();

  // --- GETTERS/SETTERS ---

  /** get the time step
   *  \return the value of t[k+1] - t[k], the current time step
   */
  const double getCurrentTimeStep() const
  {
    return h;
  };

  /** set current h (ie tk+1 - tk)
   *  \param the new value for h.
   *  Warning: if the TimeDiscretisation has been built
   *  with a complete tk vector (tdCase = 1), a call to this function
   *  will switch to tdCase = 2, ie h = newH for all future steps (until a new set of h)
   */
  void setCurrentTimeStep(double);

  /** get the value of tk at step k
   *  \return a double
   */
  inline const double getTk(unsigned int ind) const
  {
    return tk.at(ind);
  }

  /** set current h (ie tk+1 - tk)
   *  \param the new value for h.
   *  Warning: if the TimeDiscretisation has been built
   *  with a complete tk vector (tdCase = 1), a call to this function
   *  will switch to tdCase = 2, ie h = newH for all future steps (until a new set of h)
   */
  void setTk(const TkVector&);

  /** Get the current time instant value ( tk[pos] ) */
  double getCurrentTime() const
  {
    return tk[pos];
  };

  /** Get time instant value at index k+1 ( tk[pos+1] ) */
  double getNextTime() const
  {
    return tk[pos + 1];
  };

  /** get the TimeDiscretisationXML of the TimeDiscretisation
   *  \return a pointer on the TimeDiscretisationXML of the TimeDiscretisation
   */
  inline TimeDiscretisationXML* getTimeDiscretisationXMLPtr() const
  {
    return timeDiscretisationXML;
  }

  /** set the TimeDiscretisationXML of the TimeDiscretisation
   *  \param TimeDiscretisationXML* : the pointer to set the TimeDiscretisationXML
   */
  inline void setTimeDiscretisationXMLPtr(TimeDiscretisationXML* timediscrxml)
  {
    timeDiscretisationXML = timediscrxml;
  }

  /** get the Model
   *  \return a pointer to Model
   */
  inline Model* getModelPtr() const
  {
    return model;
  };

  /** Steps to next time step. */
  void increment();

  // --- OTHER FUNCTIONS ---
  /** print the discretisation data to the screen
   */
  void display() const;

  // --- XML Functions ---

  /** saves the TimeDiscretisation to the XML tree
   *  \exception RuntimeException
   */
  void saveTimeDiscretisationToXML();
};

#endif // TIMEDISCRETISATION_H
