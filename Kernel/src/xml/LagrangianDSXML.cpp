/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "LagrangianDSXML.hpp"
#include "SiconosMemoryXML.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosVector.hpp"

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

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_Q)))
    qNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_Q0)))
    q0Node = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_Q0 + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_QMEMORY)))
  {
    qMemoryNode = node;
    qMemoryXML.reset(new SiconosMemoryXML(qMemoryNode));
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_VELOCITY)))
    velocityNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_VELOCITY0)))
    velocity0Node = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_VELOCITY0 + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_VELOCITYMEMORY)))
  {
    velocityMemoryNode = node;
    velocityMemoryXML.reset(new SiconosMemoryXML(velocityMemoryNode));
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_Mass)))
    MassNode = node;
  else
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LNLDS_Mass + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_QNLINERTIA)))
    NNLNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_FINT)))
    FIntNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LNLDS_FEXT)))
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
    if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, name)))
      jacobianFIntNode[i] = node;
    name = "Jacobian";
    if (i)
      name += "QNNL";
    else
      name += "VelocityNNL";
    if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, name)))
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
    qMemoryXML->setSiconosMemoryVector(*smem.vectorMemory());
  }
  else
  {
    qMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    qMemoryXML->setSiconosMemoryVector(*smem.vectorMemory());
  }
}
void LagrangianDSXML::setVelocityMemory(const SiconosMemory& smem)
{
  if (hasVelocityMemory() == false)
  {
    velocityMemoryXML.reset(new SiconosMemoryXML(NULL, rootNode, LNLDS_VELOCITYMEMORY));
    velocityMemoryNode = velocityMemoryXML->getSiconosMemoryXMLNode();

    velocityMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    velocityMemoryXML->setSiconosMemoryVector(*smem.vectorMemory());
  }
  else
  {
    velocityMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    velocityMemoryXML->setSiconosMemoryVector(*smem.vectorMemory());
  }
}
