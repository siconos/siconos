//$Id: PeriodicBC.cpp,v 1.13 2005/03/08 12:41:36 jbarbier Exp $

#include "PeriodicBC.h"
#include "check.h"

PeriodicBC::PeriodicBC(): BoundaryCondition()
{
  this->boundaryType = PERIODICBC;
}

PeriodicBC::PeriodicBC(BoundaryConditionXML* bcxml): BoundaryCondition(bcxml)
{
  this->boundaryType = PERIODICBC;
}

PeriodicBC::~PeriodicBC()
{}

void PeriodicBC::fillBCWithBCXML()
{
  if (this->bcXML != NULL)
  {
    OUT("PeriodicBC::fillBCWithBCXML\n");
  }
  else RuntimeException::selfThrow("PeriodicBC::fillBCWithBCXML - The BoundaryConditionXML object doesn't exists");
}

void PeriodicBC::saveBCToXML()
{
  if (this->bcXML != NULL)
  {
    OUT("PeriodicBC::saveBCToXML\n");
  }
  else RuntimeException::selfThrow("PeriodicBC::saveBCToXML - The BoundaryConditionXML object doesn't exists");
}

void PeriodicBC::createBoundaryCondition(BoundaryConditionXML * bcXML)//, DynamicalSystem* ds)
{
  if (bcXML != NULL)
  {
    this->bcXML = bcXML;
    this->boundaryType = PERIODICBC;
    this->fillBCWithBCXML();
  }
  else
  {}
}


PeriodicBC* PeriodicBC::convert(BoundaryCondition* bc)
{
  cout << "PeriodicBC::convert (BoundaryCondition* bc)" << endl;
  PeriodicBC* pbc = dynamic_cast<PeriodicBC*>(bc);
  return pbc;
}

//$Log: PeriodicBC.cpp,v $
//Revision 1.13  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.12  2005/03/07 13:17:20  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.11  2005/01/31 16:26:24  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.10  2004/09/03 14:41:47  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.9  2004/08/17 15:12:43  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.8  2004/07/29 14:25:38  jbarbier
//- $Log: PeriodicBC.cpp,v $
//- Revision 1.13  2005/03/08 12:41:36  jbarbier
//- - constant variables files modified :
//- Some constants added in SiconosConst
//-
//- all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//-
//- Revision 1.12  2005/03/07 13:17:20  jbarbier
//- - new test : Ball2D, with a ball moving in a 2D system
//-
//- - another constant variables moved/refactored in XMLTagsName
//- - making uniform the name of the constant variables
//-
//- Revision 1.11  2005/01/31 16:26:24  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.10  2004/09/03 14:41:47  jbarbier
//- - new functions to create the boundary condition of the dynamical systems
//- - new functions to add an interaction to a NSDS
//- - new functions to create the relation and the non-smooth law of an interaciton
//-
//- Revision 1.9  2004/08/17 15:12:43  jbarbier
//- - methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//- createRelation and createNSLaw completed with the required attributes
//- and $Id: PeriodicBC.cpp,v 1.13 2005/03/08 12:41:36 jbarbier Exp $ added
//
