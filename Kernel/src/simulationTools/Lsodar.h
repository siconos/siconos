#ifndef Lsodar_H
#define Lsodar_H

#include "OneStepIntegrator.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "LsodarXML.h"
#include "check.h"

//#include "C_file_for_Lsodar.h"

using namespace std;

/** \class Lsodar
 *  \brief It's a kind of single-step Integrator
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 */
class Lsodar : public OneStepIntegrator
{
public:

  /** \fn Lsodar()
   *  \brief default constructor
   */
  Lsodar();

  /** \fn Lsodar(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor with XML object of the Lsodar
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   Integrator
   */
  Lsodar(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem*);

  ~Lsodar();


  /** \fn double getR()
  *   \brief Return the r of the OneStepIntegrator
  *   \return double : the value of r
  */
  inline double getR(void) const
  {
    return this->r;
  }

  /** \fn void setR(double r)
  *   \brief Return the r of OneStepIntegrator
  *   \param double : the value to set r
  */
  inline void setR(const double r)
  {
    this->r = r;
  }

  /** \fn void computeFreeState()
  *   \brief compute the free state of the dynamical system
  */
  void computeFreeState();

  /** \fn void integrate()
  *   \brief integrates the dynamical system
  */
  void integrate();

  /** \fn void saveIntegratorToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   */
  void saveIntegratorToXML();

  /** \fn void createOneStepIntegrator(OneStepIntegratorXML * osiXML, TimeDiscretisation*, DynamicalSystem*)
   *  \brief allows to create the Integrator Lsodar with an xml file, or the needed data
   *  \param OneStepNSProblemXML * : the XML object for this OneStepIntegrator
   *  \param TimeDiscretisation * : The NSDS which contains this OneStepIntegrator
   *  \exception RuntimeException
   */
  void createOneStepIntegrator(OneStepIntegratorXML * osiXML, TimeDiscretisation* td, DynamicalSystem* ds);//, Strategy * strategy = NULL);

  /** \fn void initialize()
  *  \brief initialization of the Lsodar integrator
  */
  void initialize();

  /** \fn void fillIntegratorWithIntegratorXML()
   *  \brief uses the OneStepIntegratorXML of the Lsodar Integrator to fill the fields of this Integrator
   *  \exception RuntimeException
   */
  void fillIntegratorWithIntegratorXML();

  /** \fn Lsodar* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  static Lsodar* convert(OneStepIntegrator* osi);

};
typedef void (* fctPtr)(int *sizeOfX, double *time, double *x, double *xdot);
extern "C" void tryfunction(fctPtr);

#endif // Lsodar_H
//$Log: Lsodar.h,v $
//Revision 1.10  2005/02/14 09:52:22  charlety
//_ getters / setters put inline
//
//Revision 1.9  2005/01/31 16:26:25  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.8  2004/09/22 14:11:13  charlety
//
//  _ revision of Doxygen comments in modelstrategy
//
//Revision 1.7  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.6  2004/09/10 11:26:16  charlety
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
//Revision 1.5  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.4  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.3  2004/08/10 13:04:19  charlety
//
//_ try to pass a pointer of a C function of th eplugin to a numerical routine
//
//Revision 1.2  2004/08/09 15:00:54  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.1  2004/08/05 14:33:21  charlety
//
//_ class LSODAR is now named Lsodar
//
//Revision 1.12  2004/07/29 14:25:39  jbarbier
