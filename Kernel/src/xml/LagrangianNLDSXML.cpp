
#include "LagrangianNLDSXML.h"

#include "check.h"


LagrangianNLDSXML::LagrangianNLDSXML() : DSXML()
{
  this->qMemoryXML = NULL;
  this->velocityMemoryXML = NULL;

  this->parentNode = NULL;

  this->qNode = NULL;
  this->q0Node = NULL;
  this->qMemoryNode = NULL;
  this->velocityNode = NULL;
  this->velocity0Node = NULL;
  this->velocityMemoryNode = NULL;
  this->QNLInertiaNode = NULL;
  this->FintNode = NULL;
  this->FextNode = NULL;
  this->jacobianQFintNode = NULL;
  this->jacobianVelocityFintNode = NULL;
  this->jacobianQQNLInertiaNode = NULL;
  this->jacobianVelocityQNLInertiaNode = NULL;
  this->MNode = NULL;
  this->ndofNode = NULL;
}

LagrangianNLDSXML::LagrangianNLDSXML(xmlNode * LagrangianNLDSNode, bool isBVP)
  : DSXML(LagrangianNLDSNode, isBVP)
{
  this->qMemoryXML = NULL;
  this->velocityMemoryXML = NULL;

  this->loadLagrangianNLDSProperties();
  this->parentNode = NULL;
}

LagrangianNLDSXML::~LagrangianNLDSXML()
{
  //  if( this->qMemoryXML != NULL ) delete this->qMemoryXML;
  //  if( this->velocityMemoryXML != NULL ) delete this->velocityMemoryXML;
}

void LagrangianNLDSXML::loadLagrangianNLDSProperties()
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_Q)) != NULL)
  {
    this->qNode = node;
  }
  else
  {
    cout << "LagrangianNLDSXML - loadLagrangianNLDSProperties Warning : tag " << LNLDS_Q << " not found, it must be optional." << endl;
    this->qNode = NULL;
    //XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_Q + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_Q0)) != NULL)
  {
    this->q0Node = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_Q0 + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_QMEMORY)) != NULL)
  {
    this->qMemoryNode = node;
    this->qMemoryXML = new SiconosMemoryXML(this->qMemoryNode);
  }
  else
  {
    //XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_QMEMORY + " not found.");
    cout << "LagrangianNLDSXML - loadLagrangianNLDSProperties warning : tag " << LNLDS_QMEMORY << " not found, it is an optional attribute." << endl;
    this->qMemoryNode = NULL;
    this->qMemoryXML = NULL;
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITY)) != NULL)
  {
    this->velocityNode = node;
  }
  else
  {
    cout << "LagrangianNLDSXML - loadLagrangianNLDSProperties Warning : tag " << LNLDS_VELOCITY << " not found, it must be optional." << endl;
    this->velocityNode = NULL;
    //XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_VELOCITY + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITY0)) != NULL)
  {
    this->velocity0Node = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_VELOCITY0 + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITYMEMORY)) != NULL)
  {
    this->velocityMemoryNode = node;
    this->velocityMemoryXML = new SiconosMemoryXML(this->velocityMemoryNode);
  }
  else
  {
    //XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_VELOCITYMEMORY + " not found.");
    cout << "LagrangianNLDSXML - loadLagrangianNLDSProperties warning : tag " << LNLDS_VELOCITYMEMORY << " not found, it is an optional attribute." << endl;
    this->velocityMemoryNode = NULL;
    this->velocityMemoryXML = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_QNLINERTIA)) != NULL)
  {
    this->QNLInertiaNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_QNLINERTIA + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_FINT)) != NULL)
  {
    this->FintNode = node;
  }
  else
  {
    //XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_FINT + " not found.");
    cout << "LagrangianNLDSXML - loadLagrangianNLDSProperties warning : tag " << LNLDS_FINT << " not found, it is an optional attribute." << endl;
    this->FintNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_FEXT)) != NULL)
  {
    this->FextNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_FEXT + " not found.");
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANQFINT)) != NULL)
  {
    this->jacobianQFintNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_JACOBIANQFINT + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANVELOCITYFINT)) != NULL)
  {
    this->jacobianVelocityFintNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_JACOBIANVELOCITYFINT + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANQQNLINERTIA)) != NULL)
  {
    this->jacobianQQNLInertiaNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_JACOBIANQQNLINERTIA + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA)) != NULL)
  {
    this->jacobianVelocityQNLInertiaNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_JACOBIANVELOCITYQNLINERTIA + " not found.");
  }



  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_M)) != NULL)
  {
    this->MNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_M + " not found.");
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_NDOF)) != NULL)
  {
    this->ndofNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LNLDS_NDOF + " not found.");
  }
}

void LagrangianNLDSXML::updateDynamicalSystemXML(xmlNode* rootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  IN("LagrangianNLDSXML::updateDynamicalSystemXML\n");
  this->rootDSXMLNode = rootDSXMLNode;
  this->loadDS(ds);
  OUT("LagrangianNLDSXML::updateDynamicalSystemXML\n");
}

//SiconosMemory LagrangianNLDSXML::getQMemory()
//{
//  IN("SiconosMemory LagrangianNLDSXML::getQMemory()\n");
//
//  SiconosMemoryXML smemXML(this->qMemoryNode, this->parentNode, LNLDS_QMEMORY);
//  SiconosMemory smem(&smemXML);
//
//  OUT("SiconosMemory LagrangianNLDSXML::getQMemory()\n");
//
//  return  smem;
//}
//
//SiconosMemory LagrangianNLDSXML::getVelocityMemory()
//{
//  IN("SiconosMemory LagrangianNLDSXML::getVelocityMemory()\n");
//
//  SiconosMemoryXML smemXML(this->velocityMemoryNode, this->parentNode, LNLDS_VELOCITYMEMORY);
//  SiconosMemory smem(&smemXML);
//
//  OUT("SiconosMemory LagrangianNLDSXML::getVelocityMemory()\n");
//
//  return  smem;
//}

//$Log: LagrangianNLDSXML.cpp,v $
//Revision 1.29  2005/02/02 15:54:52  jbarbier
//- sample RollingBalls added
//
//- function getArray() added to SimpleVector to return the pointer on the array of double values
//
//Revision 1.28  2005/01/24 14:33:03  jbarbier
//- OneStepNSProblem > Solver tag is available and managed in the XML part
//
//- tests added on OneStepNSProblem > Solver tag
//
//Revision 1.27  2004/09/10 08:04:50  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.26  2004/08/23 14:30:02  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.25  2004/08/09 15:00:55  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.24  2004/08/06 11:27:53  jbarbier
//- new tests with the XML and the optional attributes
//
//- tests on the save of the XML data
//
//Revision 1.23  2004/08/04 11:03:23  jbarbier
//- about the SiconosMemory : when a SiconosMemory has a maxSize greater than the
//number of steps in memory required by an integrator, the oldest SiconosVector
//are deleted
//
//- the way to initialize the SiconosMemory by the integrator has been updated to
//match with these changes
//
//Revision 1.22  2004/08/03 12:07:12  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.21  2004/07/30 14:37:15  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.20  2004/07/29 14:25:42  jbarbier
