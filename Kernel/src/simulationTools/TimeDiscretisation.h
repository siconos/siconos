/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
#ifndef TIMEDISCRETISATION_H
#define TIMEDISCRETISATION_H

#include "SimpleVector.h"
#include "TimeDiscretisationXML.h"
#include "RuntimeException.h"
#include "check.h"
#include "Strategy.h"
#include <iostream>
#include <vector>

class Strategy;

/** \class TimeDiscretisation
 *  \brief The time discretisation scheme
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.2.
 *  \date (Creation) Apr 26, 2004
 *
 * Two types of constructors:
 * - IO (at the time only xml)
 * - straightforward
 *  The construction process is always the same:
 *   - data loading (from io or as arguments of the straightforward constructor function)
 *   - available data checking using checkCase function. This sets the value of tdCase which tells us
 *     how the timeDiscretisation will be compute, depending on the input.
 *   - compute the time discretisation members, using compute(tdCase) method.
 *
 * At the time, TimeDiscretisation is linked with the strategy, but this should change (one discretisation per integrator?).
 * That is s why get/set for \f$t_0\f$ and \f$T\f$ (from model) are encapsulate.
 * setT0 function is a friend of class Model, to allow t0 changes in Model without a call to setT0 of the model, which would
 * result in an infinite loop.
 *
 * Note that we forbid any construction of a TimeDiscretisation without a given (not NULL) strategy.
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

  /** vector of time values at each step (=> size = nSteps) */
  SimpleVector* tk;

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

  /* the strategy of simulation */
  Strategy* strategy;

  /** Flags to know if pointers have been allocated inside present class constructors or not */
  bool isTkAllocatedIn;

  /** indicator of the way of construction for this time discretisation.
      Obtained with a call to checkCase(...) during construction.*/
  int tdCase;

  /** \fn TimeDiscretisation()
   *  \brief default constructor
   */
  TimeDiscretisation();

  /** \fn const int checkCase(const bool&, const bool&,const bool&,const bool&) const;
   *  \brief determine the way to compute TD values, according to which data are provided => set tdCase
   *  \param booleans indicated if tk, h, nSteps and T are present or not
   */
  void checkCase(const bool&, const bool&, const bool&, const bool&);

  /** \fn void compute(const int&);
   *  \brief to compute TD values according a case, depending on input values and provided
   *   by checkCase.
   *  \param an int to switch to the good case
   */
  void compute(const int&);

public:

  // --- CONSTRUCTORS/DESTRUCTOR ---
  // IO constructor -> xml
  /** \fn TimeDiscretisation(TimeDiscretisationXML*, Strategy*)
   *  \brief constructor with XML
   *  \param TimeDiscretisationXML* : the XML object corresponding
   *  \param Strategy* : the strategy that owns this discretisation
   */
  TimeDiscretisation(TimeDiscretisationXML*, Strategy *);

  // --- Straightforward constructors ---

  /** \fn TimeDiscretisation(SimpleVector *, Strategy*)
   *  \brief constructor with tk and strategy as given data
   *  \param pointer on  a SimpleVector that describes the discretisation
   *  \param Strategy* : the strategy that owns this discretisation
  */
  TimeDiscretisation(SimpleVector *, Strategy*);

  /** \fn TimeDiscretisation(const double& newH, const unsigned int& newNSteps, Strategy* str)
   *  \brief constructor with h, nSteps and strategy as given data
   *  \param double (h), unsigned int (nSteps)
   *  \param Strategy* : the strategy that owns this discretisation
   */
  TimeDiscretisation(const double& newH, const unsigned int& newNSteps, Strategy* str);

  /** \fn TimeDiscretisation(const unsigned int& newNSteps, Strategy* str)
   *  \brief constructor with nSteps and strategy as given data
   *  \param int (nSteps)
   *  \param Strategy* : the strategy that owns this discretisation
  */
  TimeDiscretisation(const unsigned int& newNSteps, Strategy* str);

  /** \fn TimeDiscretisation(const double& newH, Strategy* str)
   *  \brief constructor with h and strategy as given data
   *  \param double (h)
   *  \param Strategy* : the strategy that owns this discretisation
   */
  TimeDiscretisation(const double& newH, Strategy* str);

  // Destructor
  ~TimeDiscretisation();

  // --- GETTERS/SETTERS ---

  /** \fn const double getH() const
   *  \brief get the time step
   *  \return the value of h
   */
  inline const double getH() const
  {
    return h;
  };

  /** \fn void setH(const double&)
   *  \brief set the time step
   *  \param the new value for h
   */
  void setH(const double& newH);

  /** \fn const unsigned int getNSteps() const
   *  \brief get the number of time steps
   *  \return the value of nSteps
   */
  inline const unsigned int getNSteps() const
  {
    return nSteps;
  };

  /** \fn void setNSteps(const unsigned int&)
   *  \brief set the number of time steps
   *  \param the new value for nSteps
   */
  void setNSteps(const unsigned int& newNSteps);

  /** \fn  const SimpleVector getTk() const
   *  \brief get the value of tk
   *  \return SimpleVector
   */
  inline const SimpleVector getTk() const
  {
    return *tk;
  }

  /** \fn  const double getTk(const unsigned int& k) const
   *  \brief get the value of tk at step k
   *  \return a double
   */
  inline const double getTk(const unsigned int & k) const
  {
    return (*tk)(k);
  }

  /** \fn SimpleVector* getTkPtr() const
   *  \brief get tk
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getTkPtr() const
  {
    return tk;
  }

  /** \fn void setTk (const SimpleVector& newValue)
   *  \brief set the value of tk to newValue
   *  \param SimpleVector newValue
   */
  void setTk(const SimpleVector& newValue);

  /** \fn void setTkPtr(SimpleVector* newPtr)
   *  \brief set tk to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setTkPtr(SimpleVector *newPtr) ;

  /** \fn const double getHMin() const
   *  \brief get hMin
   *  \return a double
   */
  inline const double getHMin() const
  {
    return hMin;
  };

  /** \fn void setHMin(const double&)
   *  \brief set hMin
   *  \param the new value for hMin
   */
  inline void setHMin(const double& newhMin)
  {
    hMin = newhMin;
  };

  /** \fn const double getHMax() const
   *  \brief get hMax
   *  \return a double
   */
  inline const double getHMax() const
  {
    return hMax;
  };

  /** \fn void setHMax(const double&)
   *  \brief set hMax
   *  \param the new value for hMax
   */
  inline void setHMax(const double& newhMax)
  {
    hMax = newhMax;
  };

  /** \fn bool isConstant(void)
   *  \brief get the value of "constant", true if the TimeDiscretisation is constant
   *  \return a boolean
   */
  inline const bool isConstant() const
  {
    return constant;
  };

  /** \fn void setConstant(bool)
   *  \brief set the value of "constant"
   *  \param a boolean
   */
  inline void setConstant(const bool& newConstant)
  {
    constant = newConstant;
  };

  /** \fn inline TimeDiscretisationXML* getTimeDiscretisationXMLPtr()
   *  \brief get the TimeDiscretisationXML of the TimeDiscretisation
   *  \return a pointer on the TimeDiscretisationXML of the TimeDiscretisation
   */
  inline TimeDiscretisationXML* getTimeDiscretisationXMLPtr() const
  {
    return timeDiscretisationXML;
  }

  /** \fn inline void setTimeDiscretisationXMLPtr(TimeDiscretisationXML* timediscrxml)
   *  \brief set the TimeDiscretisationXML of the TimeDiscretisation
   *  \param TimeDiscretisationXML* : the pointer to set the TimeDiscretisationXML
   */
  inline void setTimeDiscretisationXMLPtr(TimeDiscretisationXML* timediscrxml)
  {
    timeDiscretisationXML = timediscrxml;
  }

  /** \fn const int getK() const
   *  \brief get the value of the current time step
   *  \return the value of k
   */
  inline const int getK() const
  {
    return k;
  }

  /** \fn void setK(const int& newValue)
   *  \brief set the value of K
   *  \param int : the new value for k
   */
  inline void setK(const int& newValue)
  {
    k = newValue;
  }

  /** \fn Strategy* getStrategyPtr(void)
   *  \brief get the strategy
   *  \return the strategy
   */
  inline Strategy* getStrategyPtr() const
  {
    return strategy;
  };

  // Getters and setters for time boundary value from model

  /** \fn const double getT0() const
   *  \brief get time min value
   *  \return the value of t0
   */
  const double getT0() const ;

  /** \fn void setT0(const double& newValue)
   *  \brief set initial time (friend function of class Model)
   *  \param double : the new value for t0
   */
  void setT0(const double& newValue);

  /** \fn const bool hasT() const
   *  \brief check if T, time max value is in the model or not
   *  \return a bool
   */
  const bool hasT() const;

  /** \fn const double getT() const
   *  \brief get time max value
   *  \return the value of T
   */
  const double getT() const;

  /** \fn void setT(const double& newValue)
   *  \brief set time max value
   *  \param double : the new value for t
   */
  void setT(const double& newValue);

  // --- OTHER FUNCTIONS ---
  /** \fn void increment()
   *  \brief time step increment
   */
  inline void increment()
  {
    k += 1;
  }

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  // --- XML Functions ---

  /** \fn void saveTimeDiscretisationToXML()
   *  \brief saves the TimeDiscretisation to the XML tree
   *  \exception RuntimeException
   */
  void saveTimeDiscretisationToXML();
};

#endif // TIMEDISCRETISATION_H
