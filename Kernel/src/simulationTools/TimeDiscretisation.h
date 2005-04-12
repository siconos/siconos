#ifndef TIMEDISCRETISATION_H
#define TIMEDISCRETISATION_H

#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "TimeDiscretisationXML.h"
#include "Strategy.h"
#include "OneStepIntegrator.h"

#include "RuntimeException.h"

class Strategy;
class OneStepIntegrator;

//using namespace std;

/** \class TimeDiscretisation
 *  \brief It represents the main data of Strategy and Integrators
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 * \todo Various constructors must be implemented for a constant TimeDiscretisation
 * giving for example one of the following triplet :
 *    -  t0,T,h --> Compute N
 *    -  t0,T,N --> Compute h
 *    -  t0,N,h --> Compute T
 * The construction from the XMLObject must follow the same rule
 *
 * \todo add the current time tk corresponding to the current step k
 *
 * \todo add a vector Memory for the previous time tk which have to be stored.
 *  The complete SiconosVector will be only used for a complete definition a priori of a variable timediscretisation.
 */
class TimeDiscretisation
{
public:

  /** \fn TimeDiscretisation()
   *  \brief default constructor
   */
  TimeDiscretisation();

  /** \fn TimeDiscretisation(TimeDiscretisationXML*, Strategy*)
   *  \brief constructor with XML object of the TimeDiscretisation
   *  \param TimeDiscretisationXML* : the XML object corresponding
   */
  TimeDiscretisation(TimeDiscretisationXML*);

  ~TimeDiscretisation();

  /** \fn Strategy* getStrategy(void)
   *  \brief get the strategy
   *  \return the strategy
   */
  inline Strategy* getStrategy(void) const
  {
    return this->strategy;
  };

  /** \fn double getH(void)
   *  \brief get the time step
   *  \return the value of h
   */
  inline double getH(void) const
  {
    return this->h;
  };


  /** \fn int getN(void)
    *  \brief get the number of time step
    *  \return the value of N
    */
  inline int getN(void) const
  {
    return this->N;
  };

  /** \fn SiconosVector* getTk(void)
    *  \brief get the time steps
    *  \return SiconosVector* : the value of tk
    */
  inline SimpleVector getTk(void) const
  {
    return this->tk;
  };

  /** \fn double getHMin(void)
    *  \brief get the min value of a time step
    *  \return the value of hMin
    */
  inline double getHMin(void) const
  {
    return this->hMin;
  };

  /** \fn double getHMax(void)
    *  \brief get the max value of a time step
    *  \return the value of hMax
    */
  inline double getHMax(void) const
  {
    return this->hMax;
  };

  /** \fn bool isConstant(void)
    *  \brief get the value of "constant" which tells if theTimeDiscretisation is constant
    *  \return the value of "constant"
    */
  inline bool isConstant(void) const
  {
    return this->constant;
  };


  /** \fn void setH(double h)
   *  \brief allow to set the value of h
   *  \param the double value to set
   */
  inline void setH(const double h)
  {
    this->h = h;
  };

  /** \fn void setN(int N)
   *  \brief allow to set the value of N
   *  \param the integer value to set
   */
  inline void setN(const int N)
  {
    this->N = N;
  };

  /** \fn void setTk(SiconosVector*)
   *  \brief allow to set the Vector tk
   *  \param SiconocVector* : the SiconosVector to set
   */
  inline void setTk(const SimpleVector& v)
  {
    this->tk = v;
  };

  /** \fn void setHMin(double hMin)
   *  \brief allow to set the value of hMin
   *  \param the double value to set
   */
  inline void setHMin(const double hMin)
  {
    this->hMin = hMin;
  };

  /** \fn void setHMin(double hMin)
   *  \brief allow to set the value of hMax
   *  \param the double value to set
   */
  inline void setHMax(const double hMax)
  {
    this->hMax = hMax;
  };

  /** \fn void setConstant(bool)
   *  \brief allow to set the value of "constant"
   *  \param a boolean value to set "constant"
   */
  inline void setConstant(const bool constant)
  {
    this->constant = constant;
  };

  /** \fn inline TimeDiscretisationXML* getTimeDiscretisationXML()
   *  \brief allows to get the TimeDiscretisationXML of the TimeDiscretisation
   *  \return a pointer on the TimeDiscretisationXML of the TimeDiscretisation
   */
  inline TimeDiscretisationXML* getTimeDiscretisationXML() const
  {
    return this->timeDiscretisationXML;
  }

  /** \fn inline void setTimeDiscretisationXML(TimeDiscretisationXML* timediscrxml)
   *  \brief allows to set the TimeDiscretisationXML of the TimeDiscretisation
   *  \param TimeDiscretisationXML* : the pointer to set the TimeDiscretisationXML
   */
  inline void setTimeDiscretisationXML(TimeDiscretisationXML* timediscrxml)
  {
    this->timeDiscretisationXML = timediscrxml;
  }

  /** \fn void setStrategy(Strategy* str)
   *  \brief allows to set the link to the Strategy of this TimeDiscretisation
   *  \param Strategy* : the pointer to set the Strategy
   */
  inline void setStrategy(Strategy* str)
  {
    this->strategy = str;
  }

  ///////////////////////////

  /** \fn int getK()
   *  \brief get the value of a k, the current step
   *  \return the value of k
   */
  inline int getK() const
  {
    return this->k;
  }

  /** \fn void setK(const int K)
   *  \brief allow to set the value of K
   *  \param int : the integer value to set
   */
  inline void setK(const int k)
  {
    this->k = k;
  }

  /** \fn void increment()
   *  \brief go to the next step. K is incremented, this the next time step
   */
  inline void increment()
  {
    this->k += 1;
  }

  ///////////////////////////

  /** \fn void init(double, double);
   *  \brief inits values of h, k, etc.
   *  \param the initial time and the final time.
   *  \exception RuntimeException
   */
  void init(double t0 = -1, double T = -1);

  /** \fn void fillTimeDiscretisationWithTimeDiscretisationXML()
   *  \brief uses the DiscretisationXML of the Discretisation to fill the fields
   *  \exception RuntimeException
   */
  void fillTimeDiscretisationWithTimeDiscretisationXML();

  /** \fn void saveTimeDiscretisationToXML()
   *  \brief saves the TimeDiscretisation to the XML tree
   *  \exception RuntimeException
   */
  void saveTimeDiscretisationToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn void createTimeDiscretisation(TimeDiscretisationXML * tdXML, double h, int N,
      SiconosVector * tk, double hMin, double hMax, bool constant, Strategy* str)
   *  \brief allows to create the TimeDiscretisation with an xml file, or the needed data
   *  \param TimeDiscretisationXML* : the XML object corresponding
   *  \param double : the value to initialize h
   *  \param int : the value to initialize N
   *  \param SiconosVector* : the value to initialize tk
   *  \param double : the value to initialize hMin
   *  \param double : the value to initialize hMax
   *  \param bool : the value to initialize constant
   *  \param Strategy* : if the TimeDiscretisation is for a OneStep problem, this is the only link to give
   *  \exception RuntimeException
   */
  void createTimeDiscretisation(TimeDiscretisationXML * tdXML, double h = -1.0, int N = -1,
                                SimpleVector *tk = NULL, double hMin = -1.0,
                                double hMax = -1.0, bool constant = true , Strategy* str = NULL); //, OneStepIntegrator* osi=NULL);

  /** \fn void checkTimeDiscretisation()
   *  \brief verifies if the triplet (t0,T,h), (t0,T,N) or (t0,h,N) is defined in the platform
   * and complete the missing element : N, h or T
   *  \exception RuntimeException
   */
  void checkTimeDiscretisation();


private:
  /** contains the time step */
  double h;

  /** contains the number of time step */
  int N;

  /** contains the time of each step */
  /*SiconosVector*/
  SimpleVector tk;

  /** contains the lowest value of a time step */
  double hMin;

  /** contains the highest value of a time step */
  double hMax;

  /** determines if all the Integrator have the same time step */
  bool constant;

  /** the XML object linked to the TimeDiscretisation to read XML data */
  TimeDiscretisationXML* timeDiscretisationXML;

  /** the current step */
  int k;

  /* the strategy of simulation */
  Strategy* strategy;
};

#endif // TIMEDISCRETISATION_H
