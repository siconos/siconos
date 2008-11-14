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
#include "FirstOrderNonLinearDSXML.h"
#include "SiconosMemoryXML.h"

using namespace std;

FirstOrderNonLinearDSXML::FirstOrderNonLinearDSXML():
  DynamicalSystemXML(), x0Node(NULL), xNode(NULL), MNode(NULL),
  fNode(NULL), jacobianXFNode(NULL), xMemoryNode(NULL)
{}

FirstOrderNonLinearDSXML::FirstOrderNonLinearDSXML(xmlNodePtr DSNode, bool isBVP):
  DynamicalSystemXML(DSNode, isBVP), x0Node(NULL), xNode(NULL), MNode(NULL),
  fNode(NULL), jacobianXFNode(NULL), xMemoryNode(NULL)
{
  xmlNodePtr node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, DS_X0)))
    x0Node = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, DS_X)))
    xNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, DS_M)))
    MNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, DS_F)))
    fNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, DS_JACOBIANXF)))
    jacobianXFNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, DS_XMEMORY)))
  {
    xMemoryNode = node;
    xMemoryXML.reset(new SiconosMemoryXML(xMemoryNode, rootNode, DS_XMEMORY));
  }
}

FirstOrderNonLinearDSXML::~FirstOrderNonLinearDSXML()
{
}

void FirstOrderNonLinearDSXML::setMMatrix(const SiconosMatrix&m)
{
  SiconosDOMTreeTools::setSiconosMatrixNodeValue(MNode, m);
}

void FirstOrderNonLinearDSXML::setMPlugin(const std::string& plugin)
{
  if (!MNode)
  {
    MNode = SiconosDOMTreeTools::createSingleNode(rootNode, DS_M);
    xmlNewProp(MNode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(MNode, "matrixPlugin", plugin);
}

void FirstOrderNonLinearDSXML::setFVector(const SiconosVector&v)
{
  SiconosDOMTreeTools::setSiconosVectorNodeValue(fNode, v);
}

void FirstOrderNonLinearDSXML::setFPlugin(const std::string& plugin)
{
  if (!fNode)
  {
    fNode = SiconosDOMTreeTools::createSingleNode(rootNode, DS_F);
    xmlNewProp(fNode, (xmlChar*)("vectorPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(fNode, "vectorPlugin", plugin);
}

void FirstOrderNonLinearDSXML::setJacobianXFMatrix(const SiconosMatrix&m)
{
  SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianXFNode, m);
}

void FirstOrderNonLinearDSXML::setJacobianXFPlugin(const std::string& plugin)
{
  if (!jacobianXFNode)
  {
    jacobianXFNode = SiconosDOMTreeTools::createSingleNode(rootNode, DS_JACOBIANXF);
    xmlNewProp(jacobianXFNode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(jacobianXFNode, "matrixPlugin", plugin);
}

void FirstOrderNonLinearDSXML::setXMemory(const SiconosMemory& smem)
{
  if (!hasXMemory())
  {
    xMemoryXML.reset(new SiconosMemoryXML(NULL, rootNode, DS_XMEMORY));
    xMemoryNode = xMemoryXML->getSiconosMemoryXMLNode();
    xMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    xMemoryXML->setSiconosMemoryVector(*smem.getVectorMemoryPtr());
  }
  else
  {
    xMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    xMemoryXML->setSiconosMemoryVector(*smem.getVectorMemoryPtr());
  }
}

