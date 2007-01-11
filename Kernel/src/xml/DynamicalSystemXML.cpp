/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
#include "DynamicalSystemXML.h"
#include "DynamicalSystem.h"
#include "SimpleMatrix.h"

using namespace std;

DynamicalSystemXML::DynamicalSystemXML():
  rootDynamicalSystemXMLNode(NULL), parentNode(NULL), x0Node(NULL), xNode(NULL),
  stepsInMemoryNode(NULL), xMemoryNode(NULL), rMemoryNode(NULL), MNode(NULL),
  fNode(NULL), jacobianXFNode(NULL), uSizeNode(NULL), uNode(NULL), TNode(NULL),
  xMemoryXML(NULL), rMemoryXML(NULL)
{}

DynamicalSystemXML::DynamicalSystemXML(xmlNodePtr DSNode, const bool& isBVP):
  rootDynamicalSystemXMLNode(DSNode), parentNode(NULL), x0Node(NULL), xNode(NULL),
  stepsInMemoryNode(NULL), xMemoryNode(NULL), rMemoryNode(NULL), MNode(NULL),
  fNode(NULL), jacobianXFNode(NULL), uSizeNode(NULL), uNode(NULL), TNode(NULL),
  xMemoryXML(NULL), rMemoryXML(NULL)
{
  xmlNodePtr node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_X0)) != NULL)
    x0Node = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_X)) != NULL)
    xNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_STEPSINMEMORY)) != NULL)
    stepsInMemoryNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_XMEMORY)) != NULL)
  {
    xMemoryNode = node;
    xMemoryXML = new SiconosMemoryXML(xMemoryNode, parentNode, DS_XMEMORY);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_RMEMORY)) != NULL)
  {
    rMemoryNode = node;
    rMemoryXML = new SiconosMemoryXML(rMemoryNode, parentNode, DS_RMEMORY);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_M)) != NULL)
    MNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_F)) != NULL)
    fNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_JACOBIANXF)) != NULL)
    jacobianXFNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, "uSize")) != NULL)
    uSizeNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_U)) != NULL)
    uNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_T)) != NULL)
    TNode = node;
}

DynamicalSystemXML::~DynamicalSystemXML()
{
  if (rMemoryXML != NULL) delete rMemoryXML;
  if (xMemoryXML != NULL) delete xMemoryXML;
}

void DynamicalSystemXML::setXMemory(const SiconosMemory& smem)
{
  if (!hasXMemory())
  {
    xMemoryXML = new SiconosMemoryXML(NULL, rootDynamicalSystemXMLNode, DS_XMEMORY);
    xMemoryNode = xMemoryXML->getSiconosMemoryXMLNode();
    xMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    xMemoryXML->setSiconosMemoryVector(smem.getVectorMemory());
  }
  else
  {
    xMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    xMemoryXML->setSiconosMemoryVector(smem.getVectorMemory());
  }
}

void DynamicalSystemXML::setRMemory(const SiconosMemory& smem)
{
  if (!hasRMemory())
  {
    rMemoryXML = new SiconosMemoryXML(NULL, rootDynamicalSystemXMLNode, DS_RMEMORY);
    rMemoryNode = rMemoryXML->getSiconosMemoryXMLNode();
    rMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    rMemoryXML->setSiconosMemoryVector(smem.getVectorMemory());
  }
  else
  {
    rMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    rMemoryXML->setSiconosMemoryVector(smem.getVectorMemory());
  }
}

void DynamicalSystemXML::setStepsInMemory(const unsigned int& nb)
{
  if (!hasStepsInMemory())
    stepsInMemoryNode = SiconosDOMTreeTools::createIntegerNode(rootDynamicalSystemXMLNode, DS_STEPSINMEMORY, nb);
  else SiconosDOMTreeTools::setIntegerContentValue(stepsInMemoryNode, nb);
}

void DynamicalSystemXML::setM(const SiconosMatrix& m)
{
  if (MNode != NULL)
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(MNode, m);
  else
    MNode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, DS_M, m);
}

void DynamicalSystemXML::setFVector(const SiconosVector&v)
{
  SiconosDOMTreeTools::setSiconosVectorNodeValue(fNode, v);
}

void DynamicalSystemXML::setFPlugin(const std::string& plugin)
{
  if (fNode == NULL)
  {
    fNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, DS_F);
    xmlNewProp(fNode, (xmlChar*)("vectorPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(fNode, "vectorPlugin", plugin);
}

void DynamicalSystemXML::setJacobianXFMatrix(const SiconosMatrix&m)
{
  SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianXFNode, m);
}

void DynamicalSystemXML::setJacobianXFPlugin(const std::string& plugin)
{
  if (jacobianXFNode == NULL)
  {
    jacobianXFNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, DS_JACOBIANXF);
    xmlNewProp(jacobianXFNode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(jacobianXFNode, "matrixPlugin", plugin);
}

void DynamicalSystemXML::setUVector(const SiconosVector& v)
{
  if (uNode != NULL)
    SiconosDOMTreeTools::setSiconosVectorNodeValue(uNode, v);
  else uNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, DS_U, v);
}

void DynamicalSystemXML::setUPlugin(const std::string& plugin)
{
  if (uNode == NULL)
  {
    uNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, "u");
    xmlNewProp(uNode, (xmlChar*)("vectorPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(uNode, "vectorPlugin", plugin);
}

void DynamicalSystemXML::setT(const SiconosMatrix &m)
{
  if (TNode != NULL)
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(TNode, m);
  else TNode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, DS_T, m);
}

void DynamicalSystemXML::setTPlugin(const std::string& plugin)
{
  if (TNode == NULL)
  {
    TNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, "T");
    xmlNewProp(TNode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(TNode, "matrixPlugin", plugin);
}
