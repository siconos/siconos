/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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

// to be deleted thanks to factories:
#include "NLinearBCXML.h"
#include "LinearBCXML.h"
#include "PeriodicBCXML.h"
#include "LinearDSIOXML.h"
#include "LagrangianDSIOXML.h"

using namespace std;

void DynamicalSystemXML::loadBoundaryConditionXML(xmlNode * rootBoundaryConditionNode)
{
  if (rootBoundaryConditionNode == NULL)  //BoundaryCondition is not defined
  {
    boundaryConditionXML = NULL;
  }
  else
  {
    xmlNodePtr node = SiconosDOMTreeTools::findNodeChild(rootBoundaryConditionNode);
    string type((char*)node->name);
    if (type == NON_LINEARBC_TAG)
      boundaryConditionXML = new NLinearBCXML(node);

    else if (type == LINEARBC_TAG)
      boundaryConditionXML = new LinearBCXML(node);

    else if (type == PERIODICBC_TAG)
      boundaryConditionXML = new PeriodicBCXML(node);

    else
      XMLException::selfThrow("DynamicalSystemXML : undefined boundary condition type : " + type);
  }
}

DynamicalSystemXML::DynamicalSystemXML():
  rootDynamicalSystemXMLNode(NULL), parentNode(NULL), nNode(NULL), x0Node(NULL), xNode(NULL),
  stepsInMemoryNode(NULL), xMemoryNode(NULL), xDotMemoryNode(NULL), rMemoryNode(NULL),
  vectorFieldNode(NULL), jacobianXNode(NULL), boundaryConditionNode(NULL), dsInputOutputNode(NULL),
  uSizeNode(NULL), uNode(NULL), TNode(NULL),
  boundaryConditionXML(NULL), xMemoryXML(NULL), xDotMemoryXML(NULL), rMemoryXML(NULL)
{}

DynamicalSystemXML::DynamicalSystemXML(xmlNodePtr DSNode, const bool& isBVP):
  rootDynamicalSystemXMLNode(DSNode), parentNode(NULL), nNode(NULL), x0Node(NULL), xNode(NULL),
  stepsInMemoryNode(NULL), xMemoryNode(NULL), xDotMemoryNode(NULL), rMemoryNode(NULL),
  vectorFieldNode(NULL), jacobianXNode(NULL), boundaryConditionNode(NULL), dsInputOutputNode(NULL),
  uSizeNode(NULL), uNode(NULL), TNode(NULL),
  boundaryConditionXML(NULL), xMemoryXML(NULL), xDotMemoryXML(NULL), rMemoryXML(NULL)
{
  xmlNodePtr node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_N)) != NULL)
    nNode = node;

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

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_XDOTMEMORY)) != NULL)
  {
    xDotMemoryNode = node;
    xDotMemoryXML = new SiconosMemoryXML(xDotMemoryNode, parentNode, DS_XDOTMEMORY);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_RMEMORY)) != NULL)
  {
    rMemoryNode = node;
    rMemoryXML = new SiconosMemoryXML(rMemoryNode, parentNode, DS_RMEMORY);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_VECTORFIELD)) != NULL)
    vectorFieldNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_JACOBIANX)) != NULL)
    jacobianXNode = node;

  // BC NOT YET IMPLEMENTED
  //   if ((node=SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, BOUNDARYCONDITION_TAG)) !=NULL)
  //     {
  //       loadBoundaryConditionXML(node);
  //       boundaryConditionNode = node;
  //     }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, "uSize")) != NULL)
    uSizeNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_U)) != NULL)
    uNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, DS_T)) != NULL)
    TNode = node;
}

DynamicalSystemXML::~DynamicalSystemXML()
{
  if (boundaryConditionXML != NULL) delete boundaryConditionXML;
  if (rMemoryXML != NULL) delete rMemoryXML;
  if (xMemoryXML != NULL) delete xMemoryXML;
  if (xDotMemoryXML != NULL) delete xDotMemoryXML;
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

void DynamicalSystemXML::setXDotMemory(const SiconosMemory& smem)
{

  if (!hasXDotMemory())
  {
    xDotMemoryXML = new SiconosMemoryXML(NULL, rootDynamicalSystemXMLNode, DS_XDOTMEMORY);
    xDotMemoryNode = xDotMemoryXML->getSiconosMemoryXMLNode();

    xDotMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    xDotMemoryXML->setSiconosMemoryVector(smem.getVectorMemory());
  }
  else
  {
    xDotMemoryXML->setSiconosMemorySize(smem.getMemorySize());
    xDotMemoryXML->setSiconosMemoryVector(smem.getVectorMemory());
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

void DynamicalSystemXML::setVectorFieldVector(const SiconosVector&v)
{
  SiconosDOMTreeTools::setSiconosVectorNodeValue(vectorFieldNode, v);
}

void DynamicalSystemXML::setVectorFieldPlugin(const std::string& plugin)
{
  if (vectorFieldNode == NULL)
  {
    vectorFieldNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, DS_VECTORFIELD);
    xmlNewProp(vectorFieldNode, (xmlChar*)("vectorPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(vectorFieldNode, "vectorPlugin", plugin);
}

void DynamicalSystemXML::setJacobianXMatrix(const SiconosMatrix&m)
{
  SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianXNode, m);
}

void DynamicalSystemXML::setJacobianXPlugin(const std::string& plugin)
{
  if (jacobianXNode == NULL)
  {
    jacobianXNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, DS_JACOBIANX);
    xmlNewProp(jacobianXNode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(jacobianXNode, "matrixPlugin", plugin);
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

void DynamicalSystemXML::updateDynamicalSystemXML(xmlNode* newRootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  rootDynamicalSystemXMLNode = newRootDSXMLNode;
  loadDynamicalSystem(ds);
}

void DynamicalSystemXML::loadDynamicalSystem(DynamicalSystem* ds)
{
  string type;
  xmlNode* node;

  if (ds->getBoundaryConditionPtr() != NULL)
  {
    type = ds->getBoundaryConditionPtr()->getType();
    node = xmlNewChild(rootDynamicalSystemXMLNode, NULL, (xmlChar*)BOUNDARYCONDITION_TAG.c_str(), NULL);
    if (type == NLINEARBC)
    {
      //xmlNewProp( node, (xmlChar*)TYPE_ATTRIBUTE.c_str(), (xmlChar*)NON_LINEARBC_TAG.c_str() );
      node = xmlNewChild(node, NULL, (xmlChar*)NON_LINEARBC_TAG.c_str(), NULL);
      boundaryConditionXML = new NLinearBCXML();

      // linkage between the DynamicalSystem and his DynamicalSystemXML
      ds->getBoundaryConditionPtr()->setBoundaryConditionXML(boundaryConditionXML);

      // creation of the DynamicalSystemXML
      static_cast<NLinearBCXML*>(boundaryConditionXML)->updateBoundaryConditionXML(node);  //, ds->getBoundaryCondition() );
    }
    else if (type == LINEARBC)
    {
      //xmlNewProp( node, (xmlChar*)TYPE_ATTRIBUTE.c_str(), (xmlChar*)LINEARBC_TAG.c_str() );
      node = xmlNewChild(node, NULL, (xmlChar*)LINEARBC_TAG.c_str(), NULL);
      boundaryConditionXML = new LinearBCXML();

      // linkage between the DynamicalSystem and his DynamicalSystemXML
      ds->getBoundaryConditionPtr()->setBoundaryConditionXML(boundaryConditionXML);

      // creation of the DynamicalSystemXML
      static_cast<LinearBCXML*>(boundaryConditionXML)->updateBoundaryConditionXML(node); //, ds->getBoundaryCondition() );
    }
    else if (type == PERIODICBC)
    {
      //xmlNewProp( node, (xmlChar*)TYPE_ATTRIBUTE.c_str(), (xmlChar*)PERIODICBC_TAG.c_str() );
      node = xmlNewChild(node, NULL, (xmlChar*)PERIODICBC_TAG.c_str(), NULL);
      boundaryConditionXML = new PeriodicBCXML();

      // linkage between the DynamicalSystem and his DynamicalSystemXML
      ds->getBoundaryConditionPtr()->setBoundaryConditionXML(boundaryConditionXML);

      // creation of the DynamicalSystemXML
      static_cast<PeriodicBCXML*>(boundaryConditionXML)->updateBoundaryConditionXML(node); //, ds->getBoundaryCondition() );
    }
    else
      XMLException::selfThrow("DynamicalSystemXML - loadDynamicalSystem error : undefined DynamicalSystem type : " + type);
  }

  if (ds->getDSInputOutputs().size() > 0)
  {
    int number;
    char num[32];
    map<int, DSInputOutputXML*>::iterator itdsio;
    xmlNodePtr dsioDefinitionNode, nsdsNode;
    DSInputOutputXML *dsioXML;

    nsdsNode = ds->getNonSmoothDynamicalSystemPtr()->getNonSmoothDynamicalSystemXMLPtr()->getNonSmoothDynamicalSystemXMLNode();
    dsioDefinitionNode = SiconosDOMTreeTools::findNodeChild((const xmlNodePtr)nsdsNode, DSINPUTOUTPUT_DEFINITION_TAG);
    if (dsioDefinitionNode == NULL)
      dsioDefinitionNode = xmlNewChild(nsdsNode, NULL, (xmlChar*)DSINPUTOUTPUT_DEFINITION_TAG.c_str(), NULL);

    for (unsigned int i = 0; i < ds->getDSInputOutputs().size(); i++)
    {
      if (ds->getDSInputOutput(i)->getDSInputOutputXML() == NULL)
      {
        number = ds->getDSInputOutput(i)->getNumber();
        sprintf(num, "%d", number);
        definedDSInputOutputNumbers.push_back(number);

        // verifies if this DSInputOutput has a number which not used
        itdsio = dsInputOutputXMLMap.find(number);
        if (itdsio == dsInputOutputXMLMap.end())
        {
          if (ds->getDSInputOutput(i)->getType() == LINEARDSIO)
          {
            node = xmlNewChild(dsioDefinitionNode, NULL, (xmlChar*)LINEAR_DSIO_TAG.c_str(), NULL);
            xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
            //xmlNewChild( node, NULL, (xmlChar*)DSINPUTOUTPUT_H.c_str(), NULL );
            dsioXML = new LinearDSIOXML();

            // linkage between the DSInputOutput and his DSInputOutputXML
            ds->getDSInputOutput(i)->setDSInputOutputXML(dsioXML);
            dsioXML->updateDSInputOutputXML(node, ds->getDSInputOutput(i));
            dsInputOutputXMLMap[number] = dsioXML;
          }
          else if (ds->getDSInputOutput(i)->getType() == NLINEARDSIO)
          {
            node = xmlNewChild(dsioDefinitionNode, NULL, (xmlChar*)NON_LINEAR_DSIO_TAG.c_str(), NULL);
            xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
            //xmlNewChild( node, NULL, (xmlChar*)DSINPUTOUTPUT_H.c_str(), NULL );
            dsioXML = new DSInputOutputXML();

            // linkage between the DSInputOutput and his DSInputOutputXML
            ds->getDSInputOutput(i)->setDSInputOutputXML(dsioXML);
            dsioXML->updateDSInputOutputXML(node, ds->getDSInputOutput(i));
            dsInputOutputXMLMap[number] = dsioXML;
          }
          else if (ds->getDSInputOutput(i)->getType() == LAGRANGIANDSIO)
          {
            node = xmlNewChild(dsioDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_DSIO_TAG.c_str(), NULL);
            xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
            //xmlNewChild( node, NULL, (xmlChar*)DSINPUTOUTPUT_H.c_str(), NULL );
            dsioXML = new LagrangianDSIOXML();

            // linkage between the DSInputOutput and his DSInputOutputXML
            ds->getDSInputOutput(i)->setDSInputOutputXML(dsioXML);
            dsioXML->updateDSInputOutputXML(node, ds->getDSInputOutput(i));
            dsInputOutputXMLMap[number] = dsioXML;
          }
          else XMLException::selfThrow("DSXML - loadDynamicalSystem | Error : the DSInputOutput type : " + ds->getDSInputOutput(i)->getType() + " doesn't exist!");

          /*  end of the save : saving the DynamicalSystem linked to this DSInputOutput */
          node = xmlNewChild(node, NULL, (xmlChar*)DS_CONCERNED.c_str(), NULL);
          node = xmlNewChild(node, NULL, (xmlChar*) DYNAMICAL_SYSTEM_TAG.c_str(), NULL);
          number = ds->getNumber();
          sprintf(num, "%d", number);
          xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
        }
        else cout << "DynamicalSystemXML - loadDynamicalSystem : the DSInputOutput type : " << ds->getDSInputOutput(i)->getType() << " already exists!" << endl;
      }
      else cout << "### DynamicalSystemXML updateDSInputOutputXML, DSIOXML != NULL" << endl;
    }
  }
}

DSInputOutputXML* DynamicalSystemXML::getDSInputOutputXML(const int& number)
{
  map<int, DSInputOutputXML*>::iterator it;

  it = dsInputOutputXMLMap.find(number);
  if (it == dsInputOutputXMLMap.end())
    return NULL;

  return dsInputOutputXMLMap[number];
}

void DynamicalSystemXML::setDSInputOutputXML(map<int, DSInputOutputXML*> m)
{
  definedDSInputOutputNumbers.clear();

  map<int, DSInputOutputXML*>::iterator iter;
  for (iter = m.begin(); iter != m.end(); iter++)
    definedDSInputOutputNumbers.push_back((*iter).first);

  dsInputOutputXMLMap = m;
}

