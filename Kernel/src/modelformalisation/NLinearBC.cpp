//$Id: NLinearBC.cpp,v 1.12 2005/03/08 12:41:36 jbarbier Exp $

#include "NLinearBC.h"
#include "check.h"

NLinearBC::NLinearBC(): BoundaryCondition()
{
  this->boundaryType = NLINEARBC;
}

NLinearBC::NLinearBC(BoundaryConditionXML* bcxml): BoundaryCondition(bcxml)
{
  this->boundaryType = NLINEARBC;
}

NLinearBC::~NLinearBC()
{}

void NLinearBC::fillBCWithBCXML()
{
  OUT("NLinearBC::fillBCWithBCXML\n");
  if (this->bcXML != NULL)
  {

  }
  else RuntimeException::selfThrow("NLinearBC::fillBCWithBCXML - The BoundaryConditionXML object doesn't exists");
}

void NLinearBC::saveBCToXML()
{
  OUT("NLinearBC::saveBCToXML\n");
  if (this->bcXML != NULL)
  {

  }
  else RuntimeException::selfThrow("NLinearBC::saveBCToXML - The BoundaryConditionXML object doesn't exists");
}

void NLinearBC::createBoundaryCondition(BoundaryConditionXML * bcXML)
{
  if (bcXML != NULL)
  {
    this->bcXML = bcXML;
    this->boundaryType = NLINEARBC;
    this->fillBCWithBCXML();
  }
  else
  {}
}

NLinearBC* NLinearBC::convert(BoundaryCondition* bc)
{
  cout << "NLinearBC::convert (BoundaryCondition* bc)" << endl;
  NLinearBC* nlbc = dynamic_cast<NLinearBC*>(bc);
  return nlbc;
}

//$Log: NLinearBC.cpp,v $
//Revision 1.12  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.11  2005/03/07 13:17:20  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.10  2005/01/31 16:26:21  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.9  2004/09/03 14:41:42  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.8  2004/08/17 15:12:40  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.7  2004/07/29 14:25:37  jbarbier
//- $Log: NLinearBC.cpp,v $
//- Revision 1.12  2005/03/08 12:41:36  jbarbier
//- - constant variables files modified :
//- Some constants added in SiconosConst
//-
//- all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//-
//- Revision 1.11  2005/03/07 13:17:20  jbarbier
//- - new test : Ball2D, with a ball moving in a 2D system
//-
//- - another constant variables moved/refactored in XMLTagsName
//- - making uniform the name of the constant variables
//-
//- Revision 1.10  2005/01/31 16:26:21  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.9  2004/09/03 14:41:42  jbarbier
//- - new functions to create the boundary condition of the dynamical systems
//- - new functions to add an interaction to a NSDS
//- - new functions to create the relation and the non-smooth law of an interaciton
//-
//- Revision 1.8  2004/08/17 15:12:40  jbarbier
//- - methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//- createRelation and createNSLaw completed with the required attributes
//- and $Id: NLinearBC.cpp,v 1.12 2005/03/08 12:41:36 jbarbier Exp $ added
//
