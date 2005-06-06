#include "LagrangianLinearTIDSXML.h"
using namespace std;

LagrangianLinearTIDSXML::LagrangianLinearTIDSXML()
  : LagrangianDSXML()
{
  this->KNode = NULL;
  this->CNode = NULL;
}

LagrangianLinearTIDSXML::LagrangianLinearTIDSXML(xmlNode * LagrangianLinearTIDSNode, bool isBVP)
// : LagrangianDSXML(LagrangianLinearTIDSNode, isBVP)
{
  this->KNode = NULL;
  this->CNode = NULL;
  this->boundaryConditionXML = NULL;
  this->rootDSXMLNode = LagrangianLinearTIDSNode;
  this->parentNode = NULL;

  this->loadLagrangianLinearTIDSProperties(isBVP);
}


void LagrangianLinearTIDSXML::loadLagrangianLinearTIDSProperties(bool isBVP)
{
  xmlNode *node;

  this->loadDSProperties(isBVP);

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_Q)) != NULL)
  {
    this->qNode = node;
  }
  else
  {
    cout << "LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties Warning : tag " << LNLDS_Q << " not found, it must be optional." << endl;
    this->qNode = NULL;
    //XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LNLDS_Q + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_Q0)) != NULL)
  {
    this->q0Node = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LNLDS_Q0 + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_QMEMORY)) != NULL)
  {
    this->qMemoryNode = node;
    this->qMemoryXML = new SiconosMemoryXML(this->qMemoryNode);
  }
  else
  {
    //XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LNLDS_QMEMORY + " not found.");
    cout << "LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties warning : tag " << LNLDS_QMEMORY << " not found, it is an optional attribute." << endl;
    this->qMemoryNode = NULL;
    this->qMemoryXML = NULL;
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITY)) != NULL)
  {
    this->velocityNode = node;
  }
  else
  {
    cout << "LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties Warning : tag " << LNLDS_VELOCITY << " not found, it must be optional." << endl;
    this->velocityNode = NULL;
    //XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LNLDS_VELOCITY + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITY0)) != NULL)
  {
    this->velocity0Node = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LNLDS_VELOCITY0 + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_VELOCITYMEMORY)) != NULL)
  {
    this->velocityMemoryNode = node;
    this->velocityMemoryXML = new SiconosMemoryXML(this->velocityMemoryNode);
  }
  else
  {
    //XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LNLDS_VELOCITYMEMORY + " not found.");
    cout << "LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties warning : tag " << LNLDS_VELOCITYMEMORY << " not found, it is an optional attribute." << endl;
    this->velocityMemoryNode = NULL;
    this->velocityMemoryXML = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_QNLINERTIA)) != NULL)
  {
    this->QNLInertiaNode = node;
  }
  else
  {
    cout << "LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties Warning : tag " << LNLDS_QNLINERTIA << " not found, it must be optional." << endl;
    this->QNLInertiaNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_FINT)) != NULL)
  {
    this->FintNode = node;
  }
  else
  {
    cout << "LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties Warning : tag " << LNLDS_FINT << " not found, it must be optional." << endl;
    this->FintNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_FEXT)) != NULL)
  {
    this->FextNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LNLDS_FEXT + " not found.");
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANQFINT)) != NULL)
  {
    this->jacobianQFintNode = node;
  }
  else
  {
    cout << "LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties Warning : tag " << LNLDS_JACOBIANQFINT << " not found, it must be optional." << endl;
    this->jacobianQFintNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANVELOCITYFINT)) != NULL)
  {
    this->jacobianVelocityFintNode = node;
  }
  else
  {
    cout << "LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties Warning : tag " << LNLDS_JACOBIANVELOCITYFINT << " not found, it must be optional." << endl;
    this->jacobianVelocityFintNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANQQNLINERTIA)) != NULL)
  {
    this->jacobianQQNLInertiaNode = node;
  }
  else
  {
    cout << "LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties Warning : tag " << LNLDS_JACOBIANQQNLINERTIA << " not found, it must be optional." << endl;
    this->jacobianQQNLInertiaNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA)) != NULL)
  {
    this->jacobianVelocityQNLInertiaNode = node;
  }
  else
  {
    cout << "LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties Warning : tag " << LNLDS_JACOBIANVELOCITYQNLINERTIA << " not found, it must be optional." << endl;
    this->jacobianVelocityQNLInertiaNode = NULL;
  }



  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_M)) != NULL)
  {
    this->MNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LNLDS_M + " not found.");
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LNLDS_NDOF)) != NULL)
  {
    this->ndofNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LNLDS_NDOF + " not found.");
  }
  /*--------------------------------------------------------------------*/
  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LTIDS_K)) != NULL)
  {
    this->KNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LTIDS_K + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LTIDS_C)) != NULL)
  {
    this->CNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianLinearTIDSXML - loadLagrangianLinearTIDSProperties error : tag " + LTIDS_C + " not found.");
  }
}

void LagrangianLinearTIDSXML::updateDynamicalSystemXML(xmlNode* rootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  IN("LagrangianLinearTIDSXML::updateDynamicalSystemXML\n");
  this->rootDSXMLNode = rootDSXMLNode;
  this->loadDS(ds);
  OUT("LagrangianLinearTIDSXML::updateDynamicalSystemXML\n");
}

