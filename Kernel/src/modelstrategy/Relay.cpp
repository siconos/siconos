//$Id: Relay.cpp,v 1.14 2005/03/08 14:23:44 jbarbier Exp $

#include "Relay.h"

Relay::Relay(): OneStepNSProblem()
{
  this->nspbType = RELAY_OSNSP;
}

Relay::Relay(OneStepNSProblemXML* osnspbxml): OneStepNSProblem(osnspbxml)
{
  this->nspbType = RELAY_OSNSP;
}

Relay::~Relay()
{}

void Relay::fillNSProblemWithNSProblemXML()
{
}

void Relay::saveNSProblemToXML()
{
  OneStepNSProblem::saveNSProblemToXML();
}

void Relay::createOneStepNSProblem(OneStepNSProblemXML * osnspbXML, Strategy * strategy)
{
  if (osnspbXML != NULL)
  {
    this->onestepnspbxml = osnspbXML;
    this->nspbType = RELAY_OSNSP;

    this->fillNSProblemWithNSProblemXML();
  }
  else
  {
    this->strategy = strategy;
    this->nspbType = RELAY_OSNSP;
    this->fillInteractionVector();
  }
  this->init();
}

Relay* Relay::convert(OneStepNSProblem* osnsp)
{
  cout << "Relay::convert (DynamicalSystem* osnsp)" << endl;
  Relay* r = dynamic_cast<Relay*>(osnsp);
  return r;
}

//$Log: Relay.cpp,v $
//Revision 1.14  2005/03/08 14:23:44  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.13  2005/01/31 16:26:26  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.12  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.11  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.10  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.9  2004/08/18 14:37:19  jbarbier
//- creation of Model, NSDS, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
//
//Revision 1.8  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.7  2004/07/29 14:25:40  jbarbier
//- $Log: Relay.cpp,v $
//- Revision 1.14  2005/03/08 14:23:44  jbarbier
//- - modification of constant variables :
//- in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//-
//- in simualtion tools, some constants have been moved to SiconosConst.h
//-
//- Revision 1.13  2005/01/31 16:26:26  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.12  2004/09/21 11:49:10  jbarbier
//- - correction in the XML save for a manual construction of the platform :
//-     DS_Concerned of the Interaction
//-     DS_Concerned of the Integrator
//-
//- - test updated for these changes
//-
//- Revision 1.11  2004/09/15 13:23:13  jbarbier
//- - corrections in the OneStepNSProblem, for the XML save. The list of interaction
//- linked to the onestepnsproblem is now saved correctly. It is updated before
//- during the creation process.
//-
//- Revision 1.10  2004/09/09 08:57:44  jbarbier
//- - functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//- createTimeDiscretisation of the Strategy done.
//-
//- => all functions to create manually the objects of the platform are done
//-
//- Revision 1.9  2004/08/18 14:37:19  jbarbier
//- - creation of Model, NSDS, Strategy(TimeStepping and EventDriven) and
//- DynamicalSystem available when the creation is in a command program
//-
//- Revision 1.8  2004/08/12 11:55:18  jbarbier
//- - new methods createModel, createNSDS, createStrategy, ...
//- they now allow to make the link with upper objects of the platform
//- it will be used for the creation of the platform without XML input file
//-
//- - the createModel method is finished but the attributes of the other objects
//- of the platform are missing for the conctruction
//- and $Id: Relay.cpp,v 1.14 2005/03/08 14:23:44 jbarbier Exp $ added
//
