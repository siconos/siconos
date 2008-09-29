/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
  DynamicalSystemXML(), qNode(NULL), q0Node(NULL), qMemoryNode(NULL), velocityNode(NULL),
  velocity0Node(NULL), velocityMemoryNode(NULL),  MassNode(NULL), NNLNode(NULL), FIntNode(NULL), FExtNode(NULL)
{}

LagrangianDSXML::LagrangianDSXML(xmlNodePtr  DSNode, bool isBVP):
  DynamicalSystemXML(DSNode, isBVP), qNode(NULL), q0Node(NULL), qMemoryNode(NULL), velocityNode(NULL),
  velocity0Node(NULL), velocityMemoryNode(NULL),  MassNode(NULL), NNLNode(NULL), FIntNode(NULL), FExtNode(NULL)
{
  xmlNodePtr node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_Q)) != NULL)
    qNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_Q0)) != NULL)
    q0Node = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_Q0 + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_QMEMORY)) != NULL)
  {
    qMemoryNode = node;
    qMemoryXML.reset(new SiconosMemoryXML(qMemoryNode));
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_VELOCITY)) != NULL)
    velocityNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_VELOCITY0)) != NULL)
    velocity0Node = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_VELOCITY0 + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_VELOCITYMEMORY)) != NULL)
  {
    velocityMemoryNode = node;
    velocityMemoryXML.reset(new SiconosMemoryXML(velocityMemoryNode));
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_Mass)) != NULL)
    MassNode = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_Mass + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_QNLINERTIA)) != NULL)
    NNLNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_FINT)) != NULL)
    FIntNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_FEXT)) != NULL)
    FExtNode = node;

  jacobianNNLNode.resize(2, NULL);
  jacobianFIntNode.resize(2, NULL);
  std::string name;
  for (unsigned int i = 0; i < 2; ++i)
  {
    name = "Jacobian";
    if (i)
      name += "QFInt";
    else
      name += "VelocityFInt";
    if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, name)) != NULL)
      jacobianFIntNode[i] = node;
    name = "Jacobian";
    if (i)
      name += "QNNL";
    else
      name += "VelocityNNL";
    if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, name)) != NULL)
      jacobianNNLNode[i] = node;
  }
}

LagrangianDSXML::~LagrangianDSXML()
{
}

void LagrangianDSXML::setQMemory(const SiconosMemory& smem)
{
  if (!hasQMemory())
  {
    qMemoryXML.reset(new SiconosMemoryXML(NULL, rootNode, LNLDS_QMEMORY));
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
    velocityMemoryXML.reset(new SiconosMemoryXML(NULL, rootNode, LNLDS_VELOCITYMEMORY));
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
