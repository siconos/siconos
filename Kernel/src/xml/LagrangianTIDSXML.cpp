//$Id: LagrangianTIDSXML.cpp,v 1.12 2005/02/02 15:54:52 jbarbier Exp $

#include "LagrangianTIDSXML.h"

LagrangianTIDSXML::LagrangianTIDSXML()
  : LagrangianNLDSXML()
{
  this->KNode = NULL;
  this->CNode = NULL;
}

LagrangianTIDSXML::LagrangianTIDSXML(xmlNode * LagrangianTIDSNode, bool isBVP)
// : LagrangianNLDSXML(LagrangianTIDSNode, isBVP)
{
  this->KNode = NULL;
  this->CNode = NULL;
  this->boundaryConditionXML = NULL;
  this->rootDSXMLNode = LagrangianTIDSNode;
  this->parentNode = NULL;

  this->loadLagrangianTIDSProperties(isBVP);
}


void LagrangianTIDSXML::loadLagrangianTIDSProperties(bool isBVP)
{
  xmlNode *node;

  this->loadDSProperties(isBVP);

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_Q)) != NULL)
  {
    this->qNode = node;
  }
  else
  {
    cout << "LagrangianTIDSXML - loadLagrangianTIDSProperties Warning : tag " << LNLDS_Q << " not found, it must be optional." << endl;
    this->qNode = NULL;
    //XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LNLDS_Q + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_Q0)) != NULL)
  {
    this->q0Node = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LNLDS_Q0 + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_QMEMORY)) != NULL)
  {
    this->qMemoryNode = node;
    this->qMemoryXML = new SiconosMemoryXML(this->qMemoryNode);
  }
  else
  {
    //XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LNLDS_QMEMORY + " not found.");
    cout << "LagrangianTIDSXML - loadLagrangianTIDSProperties warning : tag " << LNLDS_QMEMORY << " not found, it is an optional attribute." << endl;
    this->qMemoryNode = NULL;
    this->qMemoryXML = NULL;
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITY)) != NULL)
  {
    this->velocityNode = node;
  }
  else
  {
    cout << "LagrangianTIDSXML - loadLagrangianTIDSProperties Warning : tag " << LNLDS_VELOCITY << " not found, it must be optional." << endl;
    this->velocityNode = NULL;
    //XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LNLDS_VELOCITY + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITY0)) != NULL)
  {
    this->velocity0Node = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LNLDS_VELOCITY0 + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITYMEMORY)) != NULL)
  {
    this->velocityMemoryNode = node;
    this->velocityMemoryXML = new SiconosMemoryXML(this->velocityMemoryNode);
  }
  else
  {
    //XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LNLDS_VELOCITYMEMORY + " not found.");
    cout << "LagrangianTIDSXML - loadLagrangianTIDSProperties warning : tag " << LNLDS_VELOCITYMEMORY << " not found, it is an optional attribute." << endl;
    this->velocityMemoryNode = NULL;
    this->velocityMemoryXML = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_QNLINERTIA)) != NULL)
  {
    this->QNLInertiaNode = node;
  }
  else
  {
    cout << "LagrangianTIDSXML - loadLagrangianTIDSProperties Warning : tag " << LNLDS_QNLINERTIA << " not found, it must be optional." << endl;
    this->QNLInertiaNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_FINT)) != NULL)
  {
    this->FintNode = node;
  }
  else
  {
    cout << "LagrangianTIDSXML - loadLagrangianTIDSProperties Warning : tag " << LNLDS_FINT << " not found, it must be optional." << endl;
    this->FintNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_FEXT)) != NULL)
  {
    this->FextNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LNLDS_FEXT + " not found.");
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANQFINT)) != NULL)
  {
    this->jacobianQFintNode = node;
  }
  else
  {
    cout << "LagrangianTIDSXML - loadLagrangianTIDSProperties Warning : tag " << LNLDS_JACOBIANQFINT << " not found, it must be optional." << endl;
    this->jacobianQFintNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANVELOCITYFINT)) != NULL)
  {
    this->jacobianVelocityFintNode = node;
  }
  else
  {
    cout << "LagrangianTIDSXML - loadLagrangianTIDSProperties Warning : tag " << LNLDS_JACOBIANVELOCITYFINT << " not found, it must be optional." << endl;
    this->jacobianVelocityFintNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANQQNLINERTIA)) != NULL)
  {
    this->jacobianQQNLInertiaNode = node;
  }
  else
  {
    cout << "LagrangianTIDSXML - loadLagrangianTIDSProperties Warning : tag " << LNLDS_JACOBIANQQNLINERTIA << " not found, it must be optional." << endl;
    this->jacobianQQNLInertiaNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA)) != NULL)
  {
    this->jacobianVelocityQNLInertiaNode = node;
  }
  else
  {
    cout << "LagrangianTIDSXML - loadLagrangianTIDSProperties Warning : tag " << LNLDS_JACOBIANVELOCITYQNLINERTIA << " not found, it must be optional." << endl;
    this->jacobianVelocityQNLInertiaNode = NULL;
  }



  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_M)) != NULL)
  {
    this->MNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LNLDS_M + " not found.");
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_NDOF)) != NULL)
  {
    this->ndofNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LNLDS_NDOF + " not found.");
  }
  /*--------------------------------------------------------------------*/
  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LTIDS_K)) != NULL)
  {
    this->KNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LTIDS_K + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LTIDS_C)) != NULL)
  {
    this->CNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianTIDSXML - loadLagrangianTIDSProperties error : tag " + LTIDS_C + " not found.");
  }
}

void LagrangianTIDSXML::updateDynamicalSystemXML(xmlNode* rootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  IN("LagrangianTIDSXML::updateDynamicalSystemXML\n");
  this->rootDSXMLNode = rootDSXMLNode;
  this->loadDS(ds);
  OUT("LagrangianTIDSXML::updateDynamicalSystemXML\n");
}

//$Log: LagrangianTIDSXML.cpp,v $
//Revision 1.12  2005/02/02 15:54:52  jbarbier
//- sample RollingBalls added
//
//- function getArray() added to SimpleVector to return the pointer on the array of double values
//
//Revision 1.11  2004/09/10 08:04:50  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.10  2004/08/23 14:30:02  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.9  2004/07/30 14:37:15  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.8  2004/07/29 14:25:43  jbarbier
//- $Log: LagrangianTIDSXML.cpp,v $
//- Revision 1.12  2005/02/02 15:54:52  jbarbier
//- - sample RollingBalls added
//-
//- - function getArray() added to SimpleVector to return the pointer on the array of double values
//-
//- Revision 1.11  2004/09/10 08:04:50  jbarbier
//- - XML save available for BoundaryCondition and Interaction
//-
//- Revision 1.10  2004/08/23 14:30:02  jbarbier
//- - All the dynamical systems can be created in a comand program and added to a
//- NSDS. The save is OK, but the creation of the boundary conditions is not yet
//- finished.
//-
//- Revision 1.9  2004/07/30 14:37:15  jbarbier
//- - saving methods for DynamicalSystemXML and LagrangianNLDSXML
//- and $Id: LagrangianTIDSXML.cpp,v 1.12 2005/02/02 15:54:52 jbarbier Exp $ added
//
