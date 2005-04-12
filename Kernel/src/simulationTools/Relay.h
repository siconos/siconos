#ifndef RELAY_H
#define RELAY_H

#include "OneStepNSProblem.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"

//using namespace std;

/** \class Relay
 *  \brief It's a way to solve NSDS. It's used in electronics
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */
class Relay : public OneStepNSProblem
{
public:
  /** \fn Relay()
   *  \brief default constructor
   */
  Relay();

  /** \fn Relay(OneStepNSProblemXML*)
   *  \brief constructor with XML object of the Relay
   *  \param OneStepNSProblemXML* : the XML object corresponding
   */
  Relay(OneStepNSProblemXML*);

  ~Relay();

  /** \fn void fillNSProblemWithNSProblemXML()
   *  \brief uses the OneStepNSProblemXML of the OneStepNSProblem to fill
   the fields of this OneStepNSProblem
   */
  void fillNSProblemWithNSProblemXML();

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   */
  void saveNSProblemToXML();

  /** \fn void createOneStepNSProblem( OneStepNSProblemXML * osiXML, Strategy * strategy )
   *  \brief allows to create the OneStepNSProblem LCP with an xml file, or the needed data
   *  \param OneStepNSProblemXML * : the XML object for this OneStepNSProblem
   *  \param Strategy * : The NSDS which contains this OneStepNSProblem
   *  \exception RuntimeException
   */
  void createOneStepNSProblem(OneStepNSProblemXML * osiXML, Strategy * strategy = NULL);

  /** \fn Relay* convert (OneStepNSProblem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static Relay* convert(OneStepNSProblem* osnsp);
};

#endif // RELAY_H
//$Log: Relay.h,v $
//Revision 1.15  2005/01/31 16:26:27  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.14  2004/09/22 14:11:14  charlety
//
//  _ revision of Doxygen comments in modelstrategy
//
//Revision 1.13  2004/09/10 11:26:17  charlety
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
//Revision 1.12  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.11  2004/08/12 11:55:19  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.10  2004/07/29 14:25:40  jbarbier
