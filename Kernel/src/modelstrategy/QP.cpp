//$Id: QP.cpp,v 1.31 2005/03/08 14:23:44 jbarbier Exp $

#include "QP.h"
#include "QPXML.h"
#include "check.h"

QP::QP(): OneStepNSProblem()
{
  this->nspbType = QP_OSNSP;
}

QP::QP(OneStepNSProblemXML* osnspbxml): OneStepNSProblem(osnspbxml)
{
  this->nspbType = QP_OSNSP;
}


QP::~QP()
{}


SiconosMatrix* QP::getQPtr(void)
{
  return &this->Q;
}

SimpleVector* QP::getPPtr(void)
{
  return &this->p;
}



void QP::formalise(double time)
{
  //OUT("QP::formaliseOneStepNSProblem\n");
}


void QP::compute(void)
{
  //OUT("QP::computeOneStepNSProblem\n");
}


void QP::fillNSProblemWithNSProblemXML()
{
  OUT("QP::fillNSProblemWithNSProblemXML");
  OneStepNSProblem::fillNSProblemWithNSProblemXML();
  if (this->onestepnspbxml != NULL)
  {
    this->Q = (static_cast<QPXML*>(this->onestepnspbxml))->getQ();
    this->p = (static_cast<QPXML*>(this->onestepnspbxml))->getP();

    //    this->display();
  }
  else RuntimeException::selfThrow("QP::fillNSProblemWithNSProblemXML - the OneStepNSProblemXML object does not exist");
}

void QP::display() const
{
  cout << "------------------------------------------------------" << endl;
  cout << "____ data of the DynamicalSystem read from a XML file" << endl;
  cout << "| Q " << endl;
  this->Q.display();
  cout << "| p " << endl ;
  this->p.display();
  cout << "____________________________" << endl;
  cout << "------------------------------------------------------" << endl;
}

void QP::saveNSProblemToXML()
{
  OUT("QP::saveNSProblemToXML");
  OneStepNSProblem::saveNSProblemToXML();
  if (this->onestepnspbxml != NULL)
  {
    /*
     * the Q et p of the LCP can be saved by calling the saveQToXML() and savePToXML()
     */
    //(static_cast<QPXML*>(this->onestepnspbxml))->setQ( &(this->Q) );
    //(static_cast<QPXML*>(this->onestepnspbxml))->setP( &(this->p) );
  }
  else RuntimeException::selfThrow("QP::saveNSProblemToXML - the OneStepNSProblemXML object does not exist");
}

void QP::savePToXML()
{
  IN("QP::savePToXML\n");
  if (this->onestepnspbxml != NULL)
  {
    (static_cast<QPXML*>(this->onestepnspbxml))->setP(&(this->p));
  }
  else RuntimeException::selfThrow("QP::savePToXML - OneStepNSProblemXML object not exists");
  OUT("QP::savePToXML\n");
}

void QP::saveQToXML()
{
  IN("QP::saveQToXML\n");
  if (this->onestepnspbxml != NULL)
  {
    (static_cast<QPXML*>(this->onestepnspbxml))->setQ(&(this->Q));
  }
  else RuntimeException::selfThrow("QP::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("QP::saveQToXML\n");
}

void QP::createOneStepNSProblem(OneStepNSProblemXML * osnspbXML, Strategy * strategy)
{
  if (osnspbXML != NULL)
  {
    this->onestepnspbxml = osnspbXML;
    this->nspbType = QP_OSNSP;

    this->fillNSProblemWithNSProblemXML();
  }
  else
  {
    this->strategy = strategy;
    this->nspbType = QP_OSNSP;
    this->fillInteractionVector();
  }
  this->init();
}


QP* QP::convert(OneStepNSProblem* osnsp)
{
  cout << "QP::convert (DynamicalSystem* osnsp)" << endl;
  QP* qp = dynamic_cast<QP*>(osnsp);
  return qp;
}

//$Log: QP.cpp,v $
//Revision 1.31  2005/03/08 14:23:44  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.30  2005/02/14 09:52:22  charlety
//_ getters / setters put inline
//
//Revision 1.29  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.28  2005/01/31 16:26:26  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.27  2005/01/27 13:57:46  jbarbier
//- suppression of old LCP and QP structures
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
//Revision 1.24  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.23  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.22  2004/09/10 11:26:17  charlety
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
//Revision 1.21  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.20  2004/08/18 14:37:19  jbarbier
//- creation of Model, NSDS, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
//
//Revision 1.19  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.18  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.17  2004/07/29 14:25:40  jbarbier
//- $Log: QP.cpp,v $
//- Revision 1.31  2005/03/08 14:23:44  jbarbier
//- - modification of constant variables :
//- in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//-
//- in simualtion tools, some constants have been moved to SiconosConst.h
//-
//- Revision 1.30  2005/02/14 09:52:22  charlety
//- _ getters / setters put inline
//-
//- Revision 1.29  2005/02/04 14:52:44  jbarbier
//- - Rolling balls in progress (contact is detected)
//-
//- - time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//-
//- Revision 1.28  2005/01/31 16:26:26  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.27  2005/01/27 13:57:46  jbarbier
//- - suppression of old LCP and QP structures
//-
//- Revision 1.26  2004/09/23 14:09:24  jbarbier
//- - modification of the integrators, the attribute r is always optional.
//-
//- - modification of the LagrangianNonLinearR. computeInput and computeOutput are
//- required.
//-
//- Revision 1.25  2004/09/22 14:11:14  charlety
//-
//-   _ revision of Doxygen comments in modelstrategy
//-
//- Revision 1.24  2004/09/21 11:49:10  jbarbier
//- - correction in the XML save for a manual construction of the platform :
//-     DS_Concerned of the Interaction
//-     DS_Concerned of the Integrator
//-
//- - test updated for these changes
//-
//- Revision 1.23  2004/09/15 13:23:13  jbarbier
//- - corrections in the OneStepNSProblem, for the XML save. The list of interaction
//- linked to the onestepnsproblem is now saved correctly. It is updated before
//- during the creation process.
//-
//- Revision 1.22  2004/09/10 11:26:17  charlety
//-
//- _ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//-
//- _ All the tests which worked with the previous version of the vector are OK with the new version.
//-
//- _ Example SICONOS and bouncingBall are OK
//-
//- _ some comments have still to be adapted to NewSiconosVector .
//-
//- _ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//-
//- Revision 1.21  2004/09/09 08:57:44  jbarbier
//- - functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//- createTimeDiscretisation of the Strategy done.
//-
//- => all functions to create manually the objects of the platform are done
//-
//- Revision 1.20  2004/08/18 14:37:19  jbarbier
//- - creation of Model, NSDS, Strategy(TimeStepping and EventDriven) and
//- DynamicalSystem available when the creation is in a command program
//-
//- Revision 1.19  2004/08/12 11:55:18  jbarbier
//- - new methods createModel, createNSDS, createStrategy, ...
//- they now allow to make the link with upper objects of the platform
//- it will be used for the creation of the platform without XML input file
//-
//- - the createModel method is finished but the attributes of the other objects
//- of the platform are missing for the conctruction
//-
//- Revision 1.18  2004/08/11 14:43:45  jbarbier
//- - beginning of the mechanism of creation without XML input file of the objects of the platform with the
//- creatObjtect methods
//-
//- - function saveWToXML for Moreau integrator, and same specific functions to save
//- M,q and Q,p for LCP and QP
//-
//- - function to check coherency of the Model
//- and $Id: QP.cpp,v 1.31 2005/03/08 14:23:44 jbarbier Exp $ added
//
