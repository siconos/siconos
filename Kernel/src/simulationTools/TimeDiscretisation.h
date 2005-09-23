#ifndef TIMEDISCRETISATION_H
#define TIMEDISCRETISATION_H

#include "SiconosMatrix.h"
#include "SimpleVector.h"
#include "TimeDiscretisationXML.h"
#include "Strategy.h"
#include "OneStepIntegrator.h"
#include "RuntimeException.h"
#include "check.h"
#include <iostream>
#include <vector>

class Strategy;
//class OneStepIntegrator;

/** \class TimeDiscretisation
 *  \brief The time discretisation scheme
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 * Two types of constructors:
 * - IO (at the time only xml)
 * - straightforward
 * Always the same process:
 * - load the data (from io or as arguments of the straightforward constructor function)
 * - determine which case we have (function checkTimeDiscretisationCase)
 * - compute the missing data using the given ones
 * 6 possible cases:
 *   | Given data     =>  Computed data
 * 1 | \f$t_k\f$      => \f$T, nSteps, h\f$
 * 2 | \f$t_k,T\f$    => \f$nSteps, h\f$
 * 3 | \f$nSteps,h\f$ => \f$T,t_k\f$
 * 5 | \f$h,T$        => \f$nSteps,t_k\f$
 * 5 | \f$nSteps,T\f$ => \f$h,t_k\f$
 * 6 |  all =>  none
 * any other case causes an exception.
 * Note that the last case is used when recovering a saved or crash file
 * (from xml for example) that contains all the data.
 * At the time, TimeDiscretisation is linked with the strategy, but this should change (one discretisation per integrator?).
 * That is s why get/set for \f$t_0\f$ and \f$T\f$ (from model) are encapsulate.
 *
 * \todo add a vector Memory for the previous time tk which have to be stored.
 *  The complete SimpleVector will be only used for a complete definition a priori of a variable timediscretisation.
 **/

class TimeDiscretisation
{
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

  /** \fn TimeDiscretisation(const double& newH, const int& newNSteps, Strategy* str)
   *  \brief constructor with h, nSteps and strategy as given data
   *  \param double (h), int (nSteps)
   *  \param Strategy* : the strategy that owns this discretisation
   */
  TimeDiscretisation(const double& newH, const int& newNSteps, Strategy* str);

  /** \fn TimeDiscretisation(const int& newNSteps, Strategy* str)
   *  \brief constructor with nSteps and strategy as given data
   *  \param int (nSteps)
   *  \param Strategy* : the strategy that owns this discretisation
  */
  TimeDiscretisation(const int& newNSteps, Strategy* str);

  /** \fn TimeDiscretisation(const double& newH, Strategy* str)
   *  \brief constructor with h and strategy as given data
   *  \param double (h)
   *  \param Strategy* : the strategy that owns this discretisation
   */
  TimeDiscretisation(const double& newH, Strategy* str);

  // Destructor
  ~TimeDiscretisation();

  // --- Constructors related functions ---

  /** \fn const int checkTimeDiscretisationCase(const bool&, const bool&,const bool&,const bool&) const;
   *  \brief determine which type of time discretisation scheme we have depending on the provided data
   *  \param booleans indicated if tk, h, nSteps and T are present or not
   *  \return an int indicating the resulting case
   */
  const int checkTimeDiscretisationCase(const bool&, const bool&, const bool&, const bool&) const;

  /** \fn void computeTimeDiscretisation(const int&);
  *  \brief check the coherency of the scheme and compute the missing values
  *  \param an in to switch on the good case
  */
  void computeTimeDiscretisation(const int&);

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
  inline void setH(const double& newH)
  {
    h = newH;
  };

  /** \fn const int getNSteps() const
   *  \brief get the number of time steps
   *  \return the value of nSteps
   */
  inline const int getNSteps() const
  {
    return nSteps;
  };

  /** \fn void setNSteps(const int&)
   *  \brief set the number of time steps
   *  \param the new value for nSteps
   */
  inline void setNSteps(const int& newNSteps)
  {
    nSteps = newNSteps;
  };

  /** \fn  const SimpleVector getTk() const
   *  \brief get the value of tk
   *  \return SimpleVector
   */
  inline const SimpleVector getTk() const
  {
    return *tk;
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
  inline void setTk(const SimpleVector& newValue)
  {
    *tk = newValue;
  }

  /** \fn void setTkPtr(SimpleVector* newPtr)
   *  \brief set tk to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setTkPtr(SimpleVector *newPtr)
  {
    if (isTkAllocatedIn = true) delete tk;
    tk = newPtr;
    isTkAllocatedIn = false;
  }

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

  /** \fn void setStrategyPtr(Strategy* str)
   *  \brief set the link to the Strategy
   *  \param Strategy* : a pointer on Strategy
   */
  inline void setStrategyPtr(Strategy* str)
  {
    strategy = str;
  }

  // Getters and setters for time boundary value from model

  /** \fn const double getT0() const
   *  \brief get time min value
   *  \return the value of t0
   */
  const double getT0() const ;

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

private:
  /** \fn TimeDiscretisation()
   *  \brief default constructor
   */
  TimeDiscretisation();

  /** contains the time step */
  double h;

  /** contains the number of time step */
  int nSteps;

  /** contains the time of each step */
  SimpleVector* tk;

  /** contains the lowest value of a time step */
  double hMin;

  /** contains the highest value of a time step */
  double hMax;

  /** true if all the Integrator have the same time step */
  bool constant;

  /** the XML object linked to the TimeDiscretisation to read XML data */
  TimeDiscretisationXML* timeDiscretisationXML;

  /** the current step */
  int k;

  /* the strategy of simulation */
  Strategy* strategy;

  /** Flags to know if pointers have been allocated inside constructors or not */
  bool isTkAllocatedIn;
};

#endif // TIMEDISCRETISATION_H
