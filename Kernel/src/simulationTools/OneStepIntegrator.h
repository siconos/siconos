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

  /** \fn OneStepIntegrator(OneStepIntegratorXML*)
   *  \brief constructor from xml file
   *  \param OneStepIntegratorXML* : the XML object corresponding
   */
  OneStepIntegrator(OneStepIntegratorXML*);

  /** \fn OneStepIntegrator(TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor from a minimum set of data
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   */
  OneStepIntegrator(TimeDiscretisation*, DynamicalSystem*);

  //OneStepIntegrator(OneStepIntegratorXML*,TimeDiscretisation*, DynamicalSystem* );

  virtual ~OneStepIntegrator();

  // getter/setter

  /** \fn const int getR() const
   *  \brief get r value
   *  \return an int
   */
  inline const int getR() const
  {
    return r;
  };

  /** \fn void setR(const int&)
   *  \brief set r
   *  \param a ref on an int
   */
  inline void setR(const int& newR)
  {
    r = newR;
  };

  /** \fn TimeDiscretisation* getTimeDiscretisationPtr()
   *  \brief get the TimeDiscretisation pointer
   *  \return a pointer on the TimeDiscretisation
   */
  inline TimeDiscretisation* getTimeDiscretisationPtr() const
  {
    return this->timeDiscretisation;
  };

  /** \fn void setTimeDiscretisationPtr(TimeDiscretisation*)
   *  \brief set the TimeDiscretisation pointer
   *  \param a pointer on a TimeDiscretisation
   */
  inline void setTimeDiscretisationPtr(TimeDiscretisation* td)
  {
    this->timeDiscretisation = td;
  };

  /** \fn DynamicalSystem* getDynamicalSystemptr()
   *  \brief get the DynamicalSystem associated with the Integrator
   *  \return a DS
   */
  inline DynamicalSystem* getDynamicalSystemPtr() const
  {
    return ds;
  };

  /** \fn void setDynamicalSystemPtr(DynamicalSystem*)
   *  \brief set the DynamicalSystem of this Integrator
   *  \param a pointer on DynamicalSystem
   */
  inline void setDynamicalSystemPtr(DynamicalSystem* newDs)
  {
    ds = newDs;
  };

  /** \fn inline OneStepIntegratorXML* getOneStepIntegratorXMLPtr()
   *  \brief get the OneStepIntegratorXML of the OneStepIntegrator
   *  \return a pointer on the OneStepIntegratorXML of the OneStepIntegrator
   */
  inline OneStepIntegratorXML* getOneStepIntegratorXMLPtr() const
  {
    return integratorxml;
  }

  /** \fn inline void setOneStepIntegratorXMLPtr(OneStepIntegratorXML* integratorxml)
   *  \brief set the OneStepIntegratorXML of the OneStepIntegrator
   *  \param OneStepIntegratorXML* : the pointer to set the OneStepIntegratorXML
   */
  inline void setOneStepIntegratorXMLPtr(OneStepIntegratorXML* newIntegratorxml)
  {
    integratorxml = newIntegratorxml;
  }

  /** \fn inline const string getType() const
   *  \brief get the type of the OneStepIntegrator
   *  \return string : the type of the OneStepIntegrator
   */
  inline const string getType() const
  {
    return integratorType;
  }

  ///////////////////////////////

  /** \fn void initialize()
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

  /** \fn void saveIntegratorToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveIntegratorToXML();

protected:

  // -- Default constructor --
  /** \fn OneStepIntegrator()
   *  \brief default constructor
   */
  OneStepIntegrator();

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
