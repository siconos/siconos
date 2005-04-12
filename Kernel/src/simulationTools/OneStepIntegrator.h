#ifndef ONESTEPINTEGRATOR_H
#define ONESTEPINTEGRATOR_H

#include "TimeDiscretisation.h"
#include "DynamicalSystem.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "OneStepIntegratorXML.h"

#include "SiconosConst.h"


//using namespace std;

class DynamicalSystem;
class TimeDiscretisation;
class OneStepIntegratorXML;

/** \class OneStepIntegrator
 *  \brief It's the generic object which can integre a DynamicalSystem
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */
class OneStepIntegrator
{
public:
  /** \fn OneStepIntegrator()
   *  \brief default constructor
   */
  OneStepIntegrator();

  /** \fn OneStepIntegrator(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor with XML object of the OneStepIntegrator
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   */
  OneStepIntegrator(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem*);

  virtual ~OneStepIntegrator();

  // getter/setter

  /** \fn int getR(void)
   *  \brief allow to get r form the OneStepIntegrator
   *  \return the value of r
   */
  inline int getR(void) const
  {
    return this->r;
  };

  /** \fn TimeDiscretisation* getTimeDiscretisation(void)
   *  \brief allow to get the TimeDiscretisation
   *  \return a pointer on the TimeDiscretisation
   */
  inline TimeDiscretisation* getTimeDiscretisation(void) const
  {
    return this->timeDiscretisation;
  };

  /** \fn void setR(int)
   *  \brief allow to set r
   *  \param the value to set r
   */
  inline void setR(const int R)
  {
    this->r = r;
  };

  /** \fn void setTimeDiscretisation(TimeDiscretisation*)
   *  \brief allow to set the TimeDiscretisation of the OneStepIntegrator
   *  \param a pointer on a TimeDiscretisation
   */
  inline void setTimeDiscretisation(TimeDiscretisation* td)
  {
    this->timeDiscretisation = td;
  };

  /** \fn DynamicalSystem* getDynamicalSystem(void)
   *  \brief allow to get the DynamicalSystem of the Integrator
   *  \return the DtnamicalSystem attached to this Integrator
   */
  inline DynamicalSystem* getDynamicalSystem(void) const
  {
    return this->ds;
  };

  /** \fn void setDynamicalSystem(DynamicalSystem*)
   *  \brief allow to set the DynamicalSystem of this Integrator
   *  \param DynamicalSystem : the DynamicalSystem used to set the DynamicalSystem ds
   */
  inline void setDynamicalSystem(DynamicalSystem* ds)
  {
    this->ds = ds;
  };

  /** \fn inline OneStepIntegratorXML* getOneStepIntegratorXML()
   *  \brief allows to get the OneStepIntegratorXML of the OneStepIntegrator
   *  \return a pointer on the OneStepIntegratorXML of the OneStepIntegrator
   */
  inline OneStepIntegratorXML* getOneStepIntegratorXML()
  {
    return this->integratorxml;
  }

  /** \fn inline void setOneStepIntegratorXML(OneStepIntegratorXML* integratorxml)
   *  \brief allows to set the OneStepIntegratorXML of the OneStepIntegrator
   *  \param OneStepIntegratorXML* : the pointer to set the OneStepIntegratorXML
   */
  inline void setOneStepIntegratorXML(OneStepIntegratorXML* integratorxml)
  {
    this->integratorxml = integratorxml;
  }

  /** \fn inline string getType()
   *  \brief allows to get the type of the OneStepIntegrator
   *  \return string : the type of the OneStepIntegrator
   */
  inline string getType()
  {
    return this->integratorType;
  }

  ///////////////////////////////

  /** \fn void initialize(bool)
  *  \brief initialise the integrator
  */
  virtual void initialize();

  /** \fn void computeFreeState(bool)
  *  \brief set r = 0 and integrates the DS
  *  \exception RuntimeException if the integration method for this type of DS is not implemented
  *  \return void
  */
  virtual void computeFreeState();

  /** \fn void integrate(bool)
   *  \brief integre some DS
   *  \exception RuntimeException if the integration method for this type of DS is not implemented
   *  \return void
   */
  virtual void integrate();

  /** \fn void nextStep(void)
   *  \brief prepares the DynamicalSystem for the next time step push x and xDot in Memory
   *  \param to be defined
   *  \exception to be defined
   *  \return void
   */
  virtual void nextStep(void);

  /** \fn void updateState()
   *  \brief update the state of the DynamicalSystem attached to this Integrator
   */
  virtual void updateState();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  virtual  void display() const;

  /** \fn void fillIntegratorWithIntegratorXML()
   *  \brief uses the OneStepIntegratorXML of the OneStepIntegrator to fill the fields of this OneStepIntegrator
   *  \exception RuntimeException
   */
  virtual void fillIntegratorWithIntegratorXML();

  /** \fn void saveIntegratorToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveIntegratorToXML();

protected:
  /** type of the Integrator */
  string integratorType;

  /** contains the DynamicalSystem for this Integrator */
  DynamicalSystem *ds;

  /** size of the memory for the integrator */
  int r;

  /** contains the data of the time discretisation */
  TimeDiscretisation *timeDiscretisation;

  /** the XML object linked to the OneStepIntegrator to read XML data */
  OneStepIntegratorXML *integratorxml;

};

#endif // ONESTEPINTEGRATOR_H
