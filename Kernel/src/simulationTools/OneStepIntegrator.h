#ifndef ONESTEPINTEGRATOR_H
#define ONESTEPINTEGRATOR_H

#include "TimeDiscretisation.h"
#include "DynamicalSystem.h"
#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "OneStepIntegratorXML.h"
#include "Strategy.h"
#include "SiconosConst.h"
#include "check.h"
#include <iostream>
#include <vector>

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

  // --- GETTERS/SETTERS ---

  /** \fn inline const string getType() const
   *  \brief get the type of the OneStepIntegrator
   *  \return string : the type of the OneStepIntegrator
   */
  inline const std::string getType() const
  {
    return integratorType;
  }

  /** \fn inline void setType(const string&)
   *  \brief set the type of the OneStepIntegrator
   *  \return string : the type of the OneStepIntegrator
   */
  inline void setType(const std::string& newType)
  {
    integratorType = newType;
  }

  /** \fn DynamicalSystem* getDynamicalSystemPtr()
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

  /** \fn const int getSizeMem() const
   *  \brief get sizeMem value
   *  \return an int
   */
  inline const int getSizeMem() const
  {
    return sizeMem;
  };

  /** \fn void setSizeMem(const int&)
   *  \brief set sizeMem
   *  \param a ref on an int
   */
  inline void setSizeMem(const int& newValue)
  {
    sizeMem = newValue;
  };

  /** \fn TimeDiscretisation* getTimeDiscretisationPtr()
   *  \brief get the TimeDiscretisation pointer
   *  \return a pointer on the TimeDiscretisation
   */
  inline TimeDiscretisation* getTimeDiscretisationPtr() const
  {
    return timeDiscretisation;
  };

  /** \fn void setTimeDiscretisationPtr(TimeDiscretisation*)
   *  \brief set the TimeDiscretisation pointer
   *  \param a pointer on a TimeDiscretisation
   */
  inline void setTimeDiscretisationPtr(TimeDiscretisation* td)
  {
    timeDiscretisation = td;
  };

  /** \fn inline OneStepIntegratorXML* getOneStepIntegratorXMLPtr()
   *  \brief get the OneStepIntegratorXML of the OneStepIntegrator
   *  \return a pointer on the OneStepIntegratorXML of the OneStepIntegrator
   */
  inline OneStepIntegratorXML* getOneStepIntegratorXMLPtr() const
  {
    return integratorXml;
  }

  /** \fn inline void setOneStepIntegratorXMLPtr(OneStepIntegratorXML* integratorXml)
   *  \brief set the OneStepIntegratorXML of the OneStepIntegrator
   *  \param OneStepIntegratorXML* : the pointer to set the OneStepIntegratorXML
   */
  inline void setOneStepIntegratorXMLPtr(OneStepIntegratorXML* newIntegratorXml)
  {
    integratorXml = newIntegratorXml;
  }


  // --- OTHERS ... ---

  /** \fn void initialize()
   *  \brief initialise the integrator
   */
  virtual void initialize();

  /** \fn void nextStep()
    *  \brief prepares the DynamicalSystem for the next time step, push x and xDot in Memory
    *  \param to be defined
    *  \exception to be defined
    *  \return void
    */
  virtual void nextStep();

  /** \fn void computeFreeState()
   *  \brief set r = 0 and integrates the DS
   *  \exception RuntimeException if the integration method for this type of DS is not implemented
   *  \return void
   */
  virtual void computeFreeState();

  /** \fn void integrate()
   *  \brief not implemented
   */
  virtual void integrate();

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

  // -- Members --
  /** type of the Integrator */
  std::string integratorType;

  /** contains the DynamicalSystem for this Integrator */
  DynamicalSystem *ds;

  /** size of the memory for the integrator */
  int sizeMem;

  /** contains the data of the time discretisation */
  TimeDiscretisation *timeDiscretisation;

  /** the XML object linked to the OneStepIntegrator to read XML data */
  OneStepIntegratorXML *integratorXml;

};

#endif // ONESTEPINTEGRATOR_H
