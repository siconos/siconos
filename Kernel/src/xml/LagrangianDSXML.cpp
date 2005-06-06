#include "LagrangianDSXML.h"
using namespace std;

LagrangianDSXML::LagrangianDSXML() : DSXML()
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

LagrangianDSXML::LagrangianDSXML(xmlNode * LagrangianDSNode, bool isBVP)
  : DSXML(LagrangianDSNode, isBVP)
{
  this->qMemoryXML = NULL;
  this->velocityMemoryXML = NULL;

  this->loadLagrangianDSProperties();
  this->parentNode = NULL;
}

LagrangianDSXML::~LagrangianDSXML()
{
  //  if( this->qMemoryXML != NULL ) delete this->qMemoryXML;
  //  if( this->velocityMemoryXML != NULL ) delete this->velocityMemoryXML;
}

void LagrangianDSXML::loadLagrangianDSProperties()
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_Q)) != NULL)
  {
    this->qNode = node;
  }
  else
  {
    cout << "LagrangianDSXML - loadLagrangianDSProperties Warning : tag " << LNLDS_Q << " not found, it must be optional." << endl;
    this->qNode = NULL;
    //XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_Q + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_Q0)) != NULL)
  {
    this->q0Node = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_Q0 + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_QMEMORY)) != NULL)
  {
    this->qMemoryNode = node;
    this->qMemoryXML = new SiconosMemoryXML(this->qMemoryNode);
  }
  else
  {
    //XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_QMEMORY + " not found.");
    cout << "LagrangianDSXML - loadLagrangianDSProperties warning : tag " << LNLDS_QMEMORY << " not found, it is an optional attribute." << endl;
    this->qMemoryNode = NULL;
    this->qMemoryXML = NULL;
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITY)) != NULL)
  {
    this->velocityNode = node;
  }
  else
  {
    cout << "LagrangianDSXML - loadLagrangianDSProperties Warning : tag " << LNLDS_VELOCITY << " not found, it must be optional." << endl;
    this->velocityNode = NULL;
    //XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_VELOCITY + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITY0)) != NULL)
  {
    this->velocity0Node = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_VELOCITY0 + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITYMEMORY)) != NULL)
  {
    this->velocityMemoryNode = node;
    this->velocityMemoryXML = new SiconosMemoryXML(this->velocityMemoryNode);
  }
  else
  {
    //XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_VELOCITYMEMORY + " not found.");
    cout << "LagrangianDSXML - loadLagrangianDSProperties warning : tag " << LNLDS_VELOCITYMEMORY << " not found, it is an optional attribute." << endl;
    this->velocityMemoryNode = NULL;
    this->velocityMemoryXML = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_QNLINERTIA)) != NULL)
  {
    this->QNLInertiaNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_QNLINERTIA + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_FINT)) != NULL)
  {
    this->FintNode = node;
  }
  else
  {
    //XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_FINT + " not found.");
    cout << "LagrangianDSXML - loadLagrangianDSProperties warning : tag " << LNLDS_FINT << " not found, it is an optional attribute." << endl;
    this->FintNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_FEXT)) != NULL)
  {
    this->FextNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_FEXT + " not found.");
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANQFINT)) != NULL)
  {
    this->jacobianQFintNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_JACOBIANQFINT + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANVELOCITYFINT)) != NULL)
  {
    this->jacobianVelocityFintNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_JACOBIANVELOCITYFINT + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANQQNLINERTIA)) != NULL)
  {
    this->jacobianQQNLInertiaNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_JACOBIANQQNLINERTIA + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA)) != NULL)
  {
    this->jacobianVelocityQNLInertiaNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_JACOBIANVELOCITYQNLINERTIA + " not found.");
  }



  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_M)) != NULL)
  {
    this->MNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_M + " not found.");
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_NDOF)) != NULL)
  {
    this->ndofNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_NDOF + " not found.");
  }
}

void LagrangianDSXML::updateDynamicalSystemXML(xmlNode* rootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  IN("LagrangianDSXML::updateDynamicalSystemXML\n");
  this->rootDSXMLNode = rootDSXMLNode;
  this->loadDS(ds);
  OUT("LagrangianDSXML::updateDynamicalSystemXML\n");
}

//SiconosMemory LagrangianDSXML::getQMemory()
//{
//  IN("SiconosMemory LagrangianDSXML::getQMemory()\n");
//
//  SiconosMemoryXML smemXML(this->qMemoryNode, this->parentNode, LNLDS_QMEMORY);
//  SiconosMemory smem(&smemXML);
//
//  OUT("SiconosMemory LagrangianDSXML::getQMemory()\n");
//
//  return  smem;
//}
//
//SiconosMemory LagrangianDSXML::getVelocityMemory()
//{
//  IN("SiconosMemory LagrangianDSXML::getVelocityMemory()\n");
//
//  SiconosMemoryXML smemXML(this->velocityMemoryNode, this->parentNode, LNLDS_VELOCITYMEMORY);
//  SiconosMemory smem(&smemXML);
//
//  OUT("SiconosMemory LagrangianDSXML::getVelocityMemory()\n");
//
//  return  smem;
//}

