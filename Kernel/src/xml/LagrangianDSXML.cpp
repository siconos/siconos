#include "LagrangianDSXML.h"
using namespace std;

LagrangianDSXML::LagrangianDSXML() :
  DSXML(), qMemoryXML(NULL), velocityMemoryXML(NULL), qNode(NULL), q0Node(NULL), qMemoryNode(NULL), velocityNode(NULL),
  velocity0Node(NULL), velocityMemoryNode(NULL), QNLInertiaNode(NULL), FintNode(NULL), FextNode(NULL), jacobianQFintNode(NULL),
  jacobianVelocityFintNode(NULL), jacobianQQNLInertiaNode(NULL), jacobianVelocityQNLInertiaNode(NULL), MNode(NULL), ndofNode(NULL)
{}

LagrangianDSXML::LagrangianDSXML(xmlNode * LagrangianDSNode, bool isBVP)
  : DSXML(LagrangianDSNode, isBVP), qMemoryXML(NULL), velocityMemoryXML(NULL), qNode(NULL), q0Node(NULL), qMemoryNode(NULL), velocityNode(NULL),
    velocity0Node(NULL), velocityMemoryNode(NULL), QNLInertiaNode(NULL), FintNode(NULL), FextNode(NULL), jacobianQFintNode(NULL),
    jacobianVelocityFintNode(NULL), jacobianQQNLInertiaNode(NULL), jacobianVelocityQNLInertiaNode(NULL), MNode(NULL), ndofNode(NULL)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_Q)) != NULL)
    qNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_Q0)) != NULL)
    q0Node = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_Q0 + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_QMEMORY)) != NULL)
  {
    qMemoryNode = node;
    qMemoryXML = new SiconosMemoryXML(qMemoryNode);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_VELOCITY)) != NULL)
    velocityNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_VELOCITY0)) != NULL)
    velocity0Node = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_VELOCITY0 + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_VELOCITYMEMORY)) != NULL)
  {
    velocityMemoryNode = node;
    velocityMemoryXML = new SiconosMemoryXML(velocityMemoryNode);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_QNLINERTIA)) != NULL)
    QNLInertiaNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_FINT)) != NULL)
    FintNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_FEXT)) != NULL)
    FextNode = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_FEXT + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_JACOBIANQFINT)) != NULL)
    jacobianQFintNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_JACOBIANVELOCITYFINT)) != NULL)
    jacobianVelocityFintNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_JACOBIANQQNLINERTIA)) != NULL)
    jacobianQQNLInertiaNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA)) != NULL)
    jacobianVelocityQNLInertiaNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_M)) != NULL)
    MNode = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_M + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LNLDS_NDOF)) != NULL)
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
    qMemoryXML = new SiconosMemoryXML(NULL, rootDSXMLNode, LNLDS_QMEMORY);
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
    velocityMemoryXML = new SiconosMemoryXML(NULL, rootDSXMLNode, LNLDS_VELOCITYMEMORY);
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
  rootDSXMLNode = newRootDSXMLNode;
  loadDS(ds);
  OUT("LagrangianDSXML::updateDynamicalSystemXML\n");
}
