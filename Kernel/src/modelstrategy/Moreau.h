//$Id: Moreau.h,v 1.30 2005/02/14 09:52:22 charlety Exp $
#ifndef MOREAU_H
#define MOREAU_H

#include "OneStepIntegrator.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "MoreauXML.h"
#include "check.h"


const int MOREAUSTEPSINMEMORY = 1;

using namespace std;

/** \class Moreau
 *  \brief It's a kind of single-step Integrator
 *  \author Jean-Michel Barbier, Jean Baptiste Charlety, Vincent ACARY
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 * $Date: 2005/02/14 09:52:22 $
 * $Revision: 1.30 $
 * $Author: charlety $
 * $Source: /CVS/Siconos/SICONOS/src/modelstrategy/Moreau.h,v $
 *
 * \todo Add the LU Factorization of W in the initialization,  and adapt the resolution in the iteration
 */
class Moreau : public OneStepIntegrator
{
public:

  /** \fn Moreau()
   *  \brief Default constructor
   */
  Moreau();

  /** \fn Moreau(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor with XML object of the Moreau
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   */
  Moreau(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem*);

  ~Moreau();

  //getter/setter
  /** \fn SiconosMatrix getW(void)
   *  \brief allows to get the SiconosMatrix W of the Moreau's Integrator
   *  \return the SiconosMatrix W
   */
  inline SiconosMatrix getW(void) const
  {
    return this->W;
  };

  /** \fn SiconosMatrix* getWPtr(void)
   *  \brief allows to get the SiconosMatrix* W of the Moreau's Integrator
   *  \return the SiconosMatrix* W
   */
  SiconosMatrix* getWPtr(void);

  /** \fn void setW(SiconosMatrix)
   *  \brief allows to set the SiconosMatrix W of the Moreau's Integrator
   *  \param the SiconosMatrix to set W
   */
  inline void setW(const SiconosMatrix& W)
  {
    this->W = W;
  };

  /** \fn double getTheta(void)
   *  \brief allows to get the double value theta of the Moreau's Integrator
   *  \return double value theta
   */
  inline double getTheta(void) const
  {
    return this->theta;
  }

  /** \fn double getTheta(void)
   *  \brief allows to set the double value theta of the Moreau's Integrator
   *  \param double value to set
   */
  inline void setTheta(const double theta)
  {
    this->theta = theta;
  }


  /** \fn void initialize()
  *  \brief initialization of the Moreau integrator; for linear time invariant systems, we compute time invariant operator (example : W)
  *  \todo LU factorization of time invariant operator (example : W)
  */
  void initialize();

  /** \fn void computeFreeState()
  *  \brief integrates the Dynamical System linked to this integrator without boring the constraints
  */
  void computeFreeState();

  /** \fn void integrate()
  *  \brief makes computations to integrate the data of a Dynamical System with the Moreau Integrator
  */
  void integrate();

  /** \fn void updateState()
  *  \brief updates the state of the Dynamical System
  */
  void updateState();

  /** \fn void saveIntegratorToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveIntegratorToXML();

  /** \fn void saveWToXML()
   *  \brief copy the matrix W of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveWToXML();


  /** \fn void fillIntegratorWithIntegratorXML()
   *  \brief uses the OneStepIntegratorXML of the Moreau Integrator to fill the fields of this Integrator
   *  \exception RuntimeException
   */
  void fillIntegratorWithIntegratorXML();

  /** \fn display()
  *  \brief Displays the data of the Moreau's integrator
  */
  void display() const;

  /** \fn void createOneStepIntegrator(OneStepIntegratorXML * osiXML,
      TimeDiscretisation* td, DynamicalSystem* ds,
      int r, double theta)
   *  \brief allows to create the Integrator Moreau with an xml file, or the needed data
   *  \param OneStepNSProblemXML * : the XML object for this OneStepIntegrator
   *  \param double :  the value for theta of this OneStepIntegrator
   *  \param TimeDiscretisation * : The NSDS which contains this OneStepIntegrator
   *  \exception RuntimeException
   */
  void createOneStepIntegrator(OneStepIntegratorXML * osiXML,
                               TimeDiscretisation* td, DynamicalSystem* ds,
                               /*int r=-1,*/ double theta = -1.0); //, Strategy * strategy = NULL);

  /** \fn Moreau* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  static Moreau* convert(OneStepIntegrator* osi);

private:
  /** a specific matrix of the Moreau Integrator */
  SiconosMatrix W;

  /** parameter of the theta method */
  double theta;
};

#endif // MOREAU_H
//$Log: Moreau.h,v $
//Revision 1.30  2005/02/14 09:52:22  charlety
//_ getters / setters put inline
//
//Revision 1.29  2005/01/31 16:26:26  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.28  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.27  2004/09/22 14:11:14  charlety
//
//  _ revision of Doxygen comments in modelstrategy
//
//Revision 1.26  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.25  2004/09/10 11:26:16  charlety
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
//Revision 1.24  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.23  2004/08/18 14:37:19  jbarbier
//- creation of Model, NSDS, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
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
//Revision 1.20  2004/07/27 14:56:05  jbarbier
//- functions createStrategy, createTimeDiscretisation and createIntegrator done
//
//Revision 1.19  2004/06/29 10:38:39  acary
//Ajout des Tag CVS ID et Log
//Ajout de la gestion pas le constructeur XML de Theta
//