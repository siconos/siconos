//$Id: OneStepIntegrator.h,v 1.28 2005/02/14 09:52:23 charlety Exp $
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


using namespace std;

class DynamicalSystem;
class TimeDiscretisation;
class OneStepIntegratorXML;

/** \class OneStepIntegrator
 *  \brief It's the generic object which can integre a DynamicalSystem
 *  \author Jean-Michel Barbier
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 * $Date: 2005/02/14 09:52:23 $
 * $Revision: 1.28 $
 * $Author: charlety $
 * $Source: /CVS/Siconos/SICONOS/src/modelstrategy/OneStepIntegrator.h,v $
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
//$Log: OneStepIntegrator.h,v $
//Revision 1.28  2005/02/14 09:52:23  charlety
//_ getters / setters put inline
//
//Revision 1.27  2005/02/10 10:35:19  jbarbier
//- new file regrouping all the const values of the model, modelingTools and numericalStrategy
//
//- new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//
//- new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//
//Revision 1.26  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.25  2004/09/22 14:11:14  charlety
//
//  _ revision of Doxygen comments in modelstrategy
//
//Revision 1.24  2004/09/14 13:49:54  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.23  2004/09/10 11:26:16  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.22  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.21  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.20  2004/08/03 12:07:11  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.19  2004/06/29 10:38:39  acary
//Ajout des Tag CVS ID et Log
//Ajout de la gestion pas le constructeur XML de Theta
//