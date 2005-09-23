#include "LagrangianDSXML.h"
using namespace std;

LagrangianDSXML::LagrangianDSXML() :
  DynamicalSystemXML(), qMemoryXML(NULL), velocityMemoryXML(NULL), qNode(NULL), q0Node(NULL), qMemoryNode(NULL), velocityNode(NULL),
  velocity0Node(NULL), velocityMemoryNode(NULL), NNLNode(NULL), FintNode(NULL), FextNode(NULL), jacobianQFintNode(NULL),
  jacobianVelocityFintNode(NULL), jacobianQNNLNode(NULL), jacobianVelocityNNLNode(NULL), MNode(NULL), ndofNode(NULL)
{}

LagrangianDSXML::LagrangianDSXML(xmlNode * LagrangianDSNode, bool isBVP)
  : DynamicalSystemXML(LagrangianDSNode, isBVP), qMemoryXML(NULL), velocityMemoryXML(NULL), qNode(NULL), q0Node(NULL), qMemoryNode(NULL), velocityNode(NULL),
    velocity0Node(NULL), velocityMemoryNode(NULL), NNLNode(NULL), FintNode(NULL), FextNode(NULL), jacobianQFintNode(NULL),
    jacobianVelocityFintNode(NULL), jacobianQNNLNode(NULL), jacobianVelocityNNLNode(NULL), MNode(NULL), ndofNode(NULL)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_Q)) != NULL)
    qNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_Q0)) != NULL)
    q0Node = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_Q0 + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_QMEMORY)) != NULL)
  {
    qMemoryNode = node;
    qMemoryXML = new SiconosMemoryXML(qMemoryNode);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_VELOCITY)) != NULL)
    velocityNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_VELOCITY0)) != NULL)
    velocity0Node = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_VELOCITY0 + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_VELOCITYMEMORY)) != NULL)
  {
    velocityMemoryNode = node;
    velocityMemoryXML = new SiconosMemoryXML(velocityMemoryNode);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_QNLINERTIA)) != NULL)
    NNLNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_FINT)) != NULL)
    FintNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_FEXT)) != NULL)
    FextNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_JACOBIANQFINT)) != NULL)
    jacobianQFintNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_JACOBIANVELOCITYFINT)) != NULL)
    jacobianVelocityFintNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_JACOBIANQQNLINERTIA)) != NULL)
    jacobianQNNLNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA)) != NULL)
    jacobianVelocityNNLNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_M)) != NULL)
    MNode = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_M + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_NDOF)) != NULL)
    ndofNode = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_NDOF + " not found.");
}

LagrangianDSXML::~LagrangianDSXML()
{}

void LagrangianDSXML::setQMemory(SiconosMemory* smem)
{
  if (!hasQMemory())
  {
    qMemoryXML = new SiconosMemoryXML(NULL, rootDynamicalSystemXMLNode, LNLDS_QMEMORY);
    qMemoryNode = qMemoryXML->getSiconosMemoryXMLNode();
    qMemoryXML->setSiconosMemorySize(smem->getMemorySize());
    qMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
  }
  else
  {
    qMemoryXML->setSiconosMemorySize(smem->getMemorySize());
    qMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
  }
}
void LagrangianDSXML::setVelocityMemory(SiconosMemory* smem)
{
  if (hasVelocityMemory() == false)
  {
    velocityMemoryXML = new SiconosMemoryXML(NULL, rootDynamicalSystemXMLNode, LNLDS_VELOCITYMEMORY);
    velocityMemoryNode = velocityMemoryXML->getSiconosMemoryXMLNode();

    velocityMemoryXML->setSiconosMemorySize(smem->getMemorySize());
    velocityMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
  }
  else
  {
    velocityMemoryXML->setSiconosMemorySize(smem->getMemorySize());
    velocityMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
  }
}

void LagrangianDSXML::updateDynamicalSystemXML(xmlNode* newRootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  IN("LagrangianDSXML::updateDynamicalSystemXML\n");
  rootDynamicalSystemXMLNode = newRootDSXMLNode;
  loadDynamicalSystem(ds);
  OUT("LagrangianDSXML::updateDynamicalSystemXML\n");
}
