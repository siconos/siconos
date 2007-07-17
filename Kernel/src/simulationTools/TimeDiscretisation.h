/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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

#include "SimpleVector.h"

class Model;
class TimeDiscretisationXML;
class SiconosVector;
class SimpleVector;

/** A time discretisation scheme
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) Apr 26, 2004
 *
 * This object is a discretization of [t0,T] defined in Model. It can belong to the Simulation, a Sensor, an Actuator ...
 *
 * Two types of constructors:
 * - IO (at the time only xml)
 * - straightforward
 *
 *  The construction process is always the same:
 *   - data loading (from io or as arguments of the straightforward constructor function)
 *   - Sets the value of tdCase which determines how the timeDiscretisation will be initialized, depending on the input.
 *  Initialization takes place during initialize of the object that owns the timeDiscretisation (simulation, sensor ...)
 *
 *
 */
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

  /* the model linked to this discretisation (used to set t0 and possibly final T) */
  Model* model;

  /** Flags to know if pointers have been allocated inside present class constructors or not */
  bool isTkAllocatedIn;

  /** indicator of the way of construction for this time discretisation.
      Obtained with a call to checkCase(...) during construction.*/
  int tdCase;

  /** Flag to check if all data are up to date (ie if it is necessary to initialize or not) */
  bool isUpToDate;

  /** default constructor
   */
  TimeDiscretisation();

  /** determine the way to compute TD values, according to which data are provided => set tdCase
  *  \param booleans indicated if tk, h and nSteps are present or not
  */
  void setCase(bool, bool, bool);

public:

  // --- CONSTRUCTORS/DESTRUCTOR ---
  // IO constructor -> xml
  /** constructor with XML
  *  \param TimeDiscretisationXML* : the XML object corresponding
  *  \param Model* : the model that owns this discretisation
  */
  TimeDiscretisation(TimeDiscretisationXML*, Model *);

  // --- Straightforward constructors ---

  /** constructor with tk and model as given data
  *  \param pointer on  a SiconosVector that describes the discretisation
  *  \param Model* : the model that owns this discretisation
  */
  TimeDiscretisation(SiconosVector *, Model*);

  /** constructor with h, nSteps and model as given data
  *  \param double (h), unsigned int (nSteps)
  *  \param Model* : the model that owns this discretisation
  */
  TimeDiscretisation(double, unsigned int, Model*);

  /** constructor with nSteps and model as given data
  *  \param int (nSteps)
  *  \param Model* : the model that owns this discretisation
  */
  TimeDiscretisation(unsigned int, Model*);

  /** constructor with h and model as given data
  *  \param double (h)
  *  \param Model* : the model that owns this discretisation
  */
  TimeDiscretisation(double, Model*);

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

  /** get the number of time steps
  *  \return the value of nSteps
  */
  inline const unsigned int getNSteps() const
  {
    return nSteps;
  };

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
  inline const double getTk(unsigned int  k) const
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
  inline void setHMin(double newhMin)
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
  inline void setHMax(double newhMax)
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
  inline void setConstant(bool newConstant)
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

  /** get the Model
  *  \return a pointer to Model
  */
  inline Model* getModelPtr() const
  {
    return model;
  };

  // Getters and setters for time boundary value from model

  /** check if T, time max value is in the model or not
   *  \return a bool
   */
  const bool hasT() const;

  /** To compute members values, depending on the input given by user at construction (ie tdCase value)
   */
  void initialize();

  // --- OTHER FUNCTIONS ---
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
