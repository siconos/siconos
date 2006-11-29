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
/*! \file
*/
#ifndef TIMEDISCRETISATION_H
#define TIMEDISCRETISATION_H

#include "TimeDiscretisationXML.h"
#include "check.h"
#include "Simulation.h"
#include <iostream>
#include <vector>

class Simulation;
class TimeDiscretisationXML;
class SiconosVector;
class SimpleVector;


/** The time discretisation scheme
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date (Creation) Apr 26, 2004
 *
 * Two types of constructors:
 * - IO (at the time only xml)
 * - straightforward
 *
 *  The construction process is always the same:
 *   - data loading (from io or as arguments of the straightforward constructor function)
 *   - available data checking using checkCase function. This sets the value of tdCase which tells us
 *     how the timeDiscretisation will be compute, depending on the input.
 *   - compute the time discretisation members, using compute(tdCase) method.
 *
 * At the time, TimeDiscretisation is linked with the simulation, but this should change (one discretisation per integrator?).
 * That is s why get/set for \f$t_0\f$ and \f$T\f$ (from model) are encapsulate.
 * setT0 function is a friend of class Model, to allow t0 changes in Model without a call to setT0 of the model, which would
 * result in an infinite loop.
 *
 * Note that we forbid any construction of a TimeDiscretisation without a given (not NULL) simulation.
 *
 * \todo add a vector Memory for the previous time tk which have to be stored.
 *  The complete SimpleVector will be only used for a complete definition a priori of a variable timediscretisation.
 **/
class TimeDiscretisation
{
private:
  /** Time step */
  double h;

  /** Number of time step */
  unsigned int nSteps;

  /** vector of time values at each step (=> size = nSteps+1) */
  SiconosVector* tk;

  /** contains the lowest possible time step value */
  double hMin;

  /** contains the highest possible time step value */
  double hMax;

  /** true if all the Integrator have the same time step */
  bool constant;

  /** the XML object linked to the TimeDiscretisation to read XML data */
  TimeDiscretisationXML* timeDiscretisationXML;

  /** Current step number*/
  int k;

  /* the simulation of simulation */
  Simulation* simulation;

  /** Flags to know if pointers have been allocated inside present class constructors or not */
  bool isTkAllocatedIn;

  /** indicator of the way of construction for this time discretisation.
      Obtained with a call to checkCase(...) during construction.*/
  int tdCase;

  /** default constructor
  */
  TimeDiscretisation();

  /** determine the way to compute TD values, according to which data are provided => set tdCase
  *  \param booleans indicated if tk, h, nSteps and T are present or not
  */
  void checkCase(const bool, const bool, const bool, const bool);

  /** to compute TD values according a case, depending on input values and provided
  *   by checkCase.
  *  \param an int to switch to the good case
  */
  void compute(const int);

public:

  // --- CONSTRUCTORS/DESTRUCTOR ---
  // IO constructor -> xml
  /** constructor with XML
  *  \param TimeDiscretisationXML* : the XML object corresponding
  *  \param Simulation* : the simulation that owns this discretisation
  */
  TimeDiscretisation(TimeDiscretisationXML*, Simulation *);

  // --- Straightforward constructors ---

  /** constructor with tk and simulation as given data
  *  \param pointer on  a SiconosVector that describes the discretisation
  *  \param Simulation* : the simulation that owns this discretisation
  */
  TimeDiscretisation(SiconosVector *, Simulation*);

  /** constructor with h, nSteps and simulation as given data
  *  \param double (h), unsigned int (nSteps)
  *  \param Simulation* : the simulation that owns this discretisation
  */
  TimeDiscretisation(const double, const unsigned int, Simulation*);

  /** constructor with nSteps and simulation as given data
  *  \param int (nSteps)
  *  \param Simulation* : the simulation that owns this discretisation
  */
  TimeDiscretisation(const unsigned int, Simulation*);

  /** constructor with h and simulation as given data
  *  \param double (h)
  *  \param Simulation* : the simulation that owns this discretisation
  */
  TimeDiscretisation(const double, Simulation*);

  // Destructor
  ~TimeDiscretisation();

  // --- GETTERS/SETTERS ---

  /** get the time step
  *  \return the value of h
  */
  inline const double getH() const
  {
    return h;
  };

  /** set the time step
  *  \param the new value for h
  */
  void setH(const double newH);

  /** get the number of time steps
  *  \return the value of nSteps
  */
  inline const unsigned int getNSteps() const
  {
    return nSteps;
  };

  /** set the number of time steps
  *  \param the new value for nSteps
  */
  void setNSteps(const unsigned int newNSteps);

  /** get the value of tk
  *  \return SimpleVector
  */
  inline const SimpleVector getTk() const
  {
    return *tk;
  }

  /** get the value of tk at step k
  *  \return a double
  */
  inline const double getTk(const unsigned int  k) const
  {
    return (*tk)(k);
  }

  /** get tk
  *  \return pointer on a SiconosVector
  */
  inline SiconosVector* getTkPtr() const
  {
    return tk;
  }

  /** set the value of tk to newValue
  *  \param SiconosVector newValue
  */
  void setTk(const SiconosVector& newValue);

  /** set tk to pointer newPtr
  *  \param SiconosVector * newPtr
  */
  void setTkPtr(SiconosVector *newPtr) ;

  /** get hMin
  *  \return a double
  */
  inline const double getHMin() const
  {
    return hMin;
  };

  /** set hMin
  *  \param the new value for hMin
  */
  inline void setHMin(const double newhMin)
  {
    hMin = newhMin;
  };

  /** get hMax
  *  \return a double
  */
  inline const double getHMax() const
  {
    return hMax;
  };

  /** set hMax
  *  \param the new value for hMax
  */
  inline void setHMax(const double newhMax)
  {
    hMax = newhMax;
  };

  /** get the value of "constant", true if the TimeDiscretisation is constant
  *  \return a boolean
  */
  inline const bool isConstant() const
  {
    return constant;
  };

  /** set the value of "constant"
  *  \param a boolean
  */
  inline void setConstant(const bool newConstant)
  {
    constant = newConstant;
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

  /** get the value of the current time step
  *  \return the value of k
  */
  inline const int getK() const
  {
    return k;
  }

  /** set the value of K
  *  \param int : the new value for k
  */
  inline void setK(const int newValue)
  {
    k = newValue;
  }

  /** get the simulation
  *  \return the simulation
  */
  inline Simulation* getSimulationPtr() const
  {
    return simulation;
  };

  // Getters and setters for time boundary value from model

  /** get time min value
  *  \return the value of t0
  */
  const double getT0() const ;

  /** set initial time (friend function of class Model)
  *  \param double : the new value for t0
  */
  void setT0(const double newValue);

  /** check if T, time max value is in the model or not
  *  \return a bool
  */
  const bool hasT() const;

  /** get time max value
  *  \return the value of T
  */
  const double getT() const;

  /** set time max value
  *  \param double : the new value for t
  */
  void setT(const double newValue);

  // --- OTHER FUNCTIONS ---
  /** time step increment
  */
  inline void increment()
  {
    k += 1;
  }

  /** print the data to the screen
  */
  void display() const;

  // --- XML Functions ---

  /** saves the TimeDiscretisation to the XML tree
  *  \exception RuntimeException
  */
  void saveTimeDiscretisationToXML();
};

#endif // TIMEDISCRETISATION_H
