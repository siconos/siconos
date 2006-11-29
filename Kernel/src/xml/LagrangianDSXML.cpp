/* Siconos-Kernel version 2.0.0, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "LagrangianDSXML.h"
#include "SiconosMemory.h"

using namespace std;

LagrangianDSXML::LagrangianDSXML() :
  DynamicalSystemXML(), ndofNode(NULL), qNode(NULL), q0Node(NULL), qMemoryNode(NULL), velocityNode(NULL),
  velocity0Node(NULL), velocityMemoryNode(NULL),  MNode(NULL), NNLNode(NULL), FintNode(NULL), FextNode(NULL), jacobianQFintNode(NULL),
  jacobianVelocityFintNode(NULL), jacobianQNNLNode(NULL), jacobianVelocityNNLNode(NULL), qMemoryXML(NULL), velocityMemoryXML(NULL)
{}

LagrangianDSXML::LagrangianDSXML(xmlNodePtr  LagrangianDSNode, const bool& isBVP)
  : DynamicalSystemXML(LagrangianDSNode, isBVP), ndofNode(NULL), qNode(NULL), q0Node(NULL), qMemoryNode(NULL), velocityNode(NULL),
    velocity0Node(NULL), velocityMemoryNode(NULL),  MNode(NULL), NNLNode(NULL), FintNode(NULL), FextNode(NULL), jacobianQFintNode(NULL),
    jacobianVelocityFintNode(NULL), jacobianQNNLNode(NULL), jacobianVelocityNNLNode(NULL), qMemoryXML(NULL), velocityMemoryXML(NULL)
{
  xmlNodePtr node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_NDOF)) != NULL)
    ndofNode = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_NDOF + " not found.");

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

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LNLDS_M)) != NULL)
    MNode = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_M + " not found.");

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
}

LagrangianDSXML::~LagrangianDSXML()
{
  if (velocityMemoryXML != NULL) delete velocityMemoryXML;
  if (qMemoryXML != NULL) delete qMemoryXML;
  qMemoryXML = NULL;
  velocityMemoryXML = NULL;
}

void LagrangianDSXML::setQMemory(const SiconosMemory& smem)
{
  if (!hasQMemory())
  {
    qMemoryXML = new SiconosMemoryXML(NULL, rootDynamicalSystemXMLNode, LNLDS_QMEMORY);
    qMemoryNode = qMemoryXML->getSiconosMemoryXMLNode();
    qMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    qMemoryXML->setSiconosMemoryVector(smem.getVectorMemory());
  }
  else
  {
    qMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    qMemoryXML->setSiconosMemoryVector(smem.getVectorMemory());
  }
}
void LagrangianDSXML::setVelocityMemory(const SiconosMemory& smem)
{
  if (hasVelocityMemory() == false)
  {
    velocityMemoryXML = new SiconosMemoryXML(NULL, rootDynamicalSystemXMLNode, LNLDS_VELOCITYMEMORY);
    velocityMemoryNode = velocityMemoryXML->getSiconosMemoryXMLNode();

    velocityMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    velocityMemoryXML->setSiconosMemoryVector(smem.getVectorMemory());
  }
  else
  {
    velocityMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    velocityMemoryXML->setSiconosMemoryVector(smem.getVectorMemory());
  }
}
