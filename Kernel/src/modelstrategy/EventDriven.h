#ifndef EVENTDRIVEN_H
#define EVENTDRIVEN_H

#include "Strategy.h"
#include "StrategyXML.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"

using namespace std;

/** \class EventDriven
 *  \brief It's a way to manage a simulation, where the event a more important than the time
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */

class EventDriven : public Strategy
{
public:
  /** \fn EventDriven()
   *  \brief defaut constructor
   */
  EventDriven();

  /** \fn EventDriven(StrategyXML*, Model*)
   *  \brief constructor with XML object of the EventDriven
   *  \param StrategyXML* : the XML object corresponding
   *  \param Model* : the Model which contains the Strategy
   */
  EventDriven(StrategyXML*, Model*);

  ~EventDriven();

  /** \fn void createStrategy(StrategyXML * strategyXML, Model * model, TimeDiscretisation * timediscretisation)
   *  \brief allows to create the Strategy with an xml file, or the needed data
   *  \param StrategyXML* : the StrategyXML linked to this Strategy
   *  \param Model* : the Model which contains this Strategy
   *  \exception RuntimeException
   */
  void createStrategy(StrategyXML * strategyXML, Model * model);//, TimeDiscretisation * timediscretisation=NULL);

  /** \fn EventDriven* convert (Strategy* str)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Strategy* : the Strategy which must be converted
   * \return a pointer on the Strategy if it is of the right type, NULL otherwise
   */
  static EventDriven* convert(Strategy* str);

private:

};

#endif // EVENTDRIVEN_H
//$Log: EventDriven.h,v $
//Revision 1.17  2005/01/31 16:26:25  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.16  2004/09/22 14:11:13  charlety
//
//  _ revision of Doxygen comments in modelstrategy
//
//Revision 1.15  2004/09/10 11:26:16  charlety
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
//Revision 1.14  2004/09/07 07:31:11  jbarbier
//- create strategies methods of the Model done
//
//Revision 1.13  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.12  2004/07/29 14:25:39  jbarbier
