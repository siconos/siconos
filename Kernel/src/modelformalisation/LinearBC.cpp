//$Id: LinearBC.cpp,v 1.17 2005/03/08 12:41:36 jbarbier Exp $

#include "LinearBC.h"
#include "check.h"

LinearBC::LinearBC(): BoundaryCondition()
{
  IN("LinearBC::LinearBC()\n");
  this->boundaryType = LINEARBC;
  this->omega = /*SiconosVector*/SimpleVector::/*SiconosVector*/SimpleVector();
  this->omega0 = SiconosMatrix::SiconosMatrix();
  this->omegaT = SiconosMatrix::SiconosMatrix();
  OUT("LinearBC::LinearBC()\n");
}

LinearBC::LinearBC(BoundaryConditionXML* bcxml): BoundaryCondition(bcxml)
{
  this->boundaryType = LINEARBC;
}

LinearBC::~LinearBC()
{}


void LinearBC::fillBCWithBCXML()
{
  OUT("LinearBC::fillBCWithBCXML\n");
  if (this->bcXML != NULL)
  {
    this->omega = static_cast<LinearBCXML*>(this->bcXML)->getOmega();
    this->omegaT = static_cast<LinearBCXML*>(this->bcXML)->getOmegaT();
    this->omega0 = static_cast<LinearBCXML*>(this->bcXML)->getOmega0();
  }
  else RuntimeException::selfThrow("LinearBC::fillBCWithBCXML - The BoundaryConditionXML object doesn't exists");
}

void LinearBC::saveBCToXML()
{
  OUT("LinearBC::saveBCToXML\n");
  if (this->bcXML != NULL)
  {
    static_cast<LinearBCXML*>(this->bcXML)->setOmega(&(this->omega));
    static_cast<LinearBCXML*>(this->bcXML)->setOmegaT(&(this->omegaT));
    static_cast<LinearBCXML*>(this->bcXML)->setOmega0(&(this->omega0));
  }
  else RuntimeException::selfThrow("LinearBC::saveBCToXML - The BoundaryConditionXML object doesn't exists");
}

void LinearBC::createBoundaryCondition(BoundaryConditionXML * bcXML,
                                       SiconosVector* omega, SiconosMatrix* omega0, SiconosMatrix* omegaT)
{
  IN("LinearBC::createBoundaryCondition\n");
  if (bcXML != NULL)
  {
    this->bcXML = bcXML;
    this->boundaryType = LINEARBC;
    this->fillBCWithBCXML();
  }
  else if (omega != NULL && omega0 != NULL && omegaT != NULL)
  {
    this->omega = *omega;
    this->omega0 = *omega0;
    this->omegaT = *omegaT;
  }
  else RuntimeException::selfThrow("LinearBC::createBoundaryCondition - The omega, omega0 and/or omegaT matrices is/are missing");
  OUT("LinearBC::createBoundaryCondition\n");
}


LinearBC* LinearBC::convert(BoundaryCondition* bc)
{
  cout << "LinearBC::convert (BoundaryCondition* bc)" << endl;
  LinearBC* lbc = dynamic_cast<LinearBC*>(bc);
  return lbc;
}

//$Log: LinearBC.cpp,v $
//Revision 1.17  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.16  2005/03/07 13:17:19  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.15  2005/02/11 17:36:01  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.14  2005/01/31 16:26:20  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.13  2005/01/18 17:07:40  charlety
//
//_ added autotools makefiles for sample directory
//
//Revision 1.12  2004/09/10 11:26:14  charlety
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
//Revision 1.11  2004/09/03 14:41:42  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.10  2004/08/17 15:12:39  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.9  2004/07/29 14:25:36  jbarbier
//- $Log: LinearBC.cpp,v $
//- Revision 1.17  2005/03/08 12:41:36  jbarbier
//- - constant variables files modified :
//- Some constants added in SiconosConst
//-
//- all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//-
//- Revision 1.16  2005/03/07 13:17:19  jbarbier
//- - new test : Ball2D, with a ball moving in a 2D system
//-
//- - another constant variables moved/refactored in XMLTagsName
//- - making uniform the name of the constant variables
//-
//- Revision 1.15  2005/02/11 17:36:01  charlety
//-
//- _ little "inspection of code"
//- _ basic getters and setters passed inline
//- _ getters functions passed const
//-
//- Revision 1.14  2005/01/31 16:26:20  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.13  2005/01/18 17:07:40  charlety
//-
//- _ added autotools makefiles for sample directory
//-
//- Revision 1.12  2004/09/10 11:26:14  charlety
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
//- Revision 1.11  2004/09/03 14:41:42  jbarbier
//- - new functions to create the boundary condition of the dynamical systems
//- - new functions to add an interaction to a NSDS
//- - new functions to create the relation and the non-smooth law of an interaciton
//-
//- Revision 1.10  2004/08/17 15:12:39  jbarbier
//- - methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//- createRelation and createNSLaw completed with the required attributes
//- and $Id: LinearBC.cpp,v 1.17 2005/03/08 12:41:36 jbarbier Exp $ added
//
