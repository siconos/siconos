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
#include "StrategyXML.h"

// includes to be deleted thanks to factories
#include "LsodarXML.h"
#include "MoreauXML.h"
#include "LCPXML.h"
#include "QPXML.h"
#include "FrictionContactXML.h"

using namespace std;

StrategyXML::StrategyXML():
  strategyNode(NULL), oneStepNSProblemXML(NULL), timeDiscretisationXML(NULL)
{}

StrategyXML::StrategyXML(xmlNode * rootStrategyNode, vector<int> definedNumberDS, vector<int> definedNumberInteraction):
  strategyNode(rootStrategyNode), oneStepNSProblemXML(NULL), timeDiscretisationXML(NULL)
{
  // Fill map of available DS
  unsigned int size = definedNumberDS.size();
  for (unsigned int i = 0; i < size; i++)
    DSAvailabilityMap[definedNumberDS[i]] = true;

  definedNumberInteractionVector = definedNumberInteraction;

  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(strategyNode, LMGC90_STRATEGY_TAG)) == NULL)
  {
    // === TimeDiscretisation data loading ===
    if ((node = SiconosDOMTreeTools::findNodeChild(strategyNode, TIMEDISCRETISATION_TAG)) != NULL)
      timeDiscretisationXML = new TimeDiscretisationXML(node);
    else
      XMLException::selfThrow("StrategyXML - strategy XML constructor  ERROR : tag " + TIMEDISCRETISATION_TAG + " not found.");

    // === OneStepIntegrator data loading ===
    if ((node = SiconosDOMTreeTools::findNodeChild(strategyNode, ONESTEPINTEGRATOR_DEFINITION_TAG)) != NULL)
    {
      xmlNode *OSInode;
      string type; // OneStepIntegrator type
      if ((OSInode = SiconosDOMTreeTools::findNodeChild((const xmlNode*)node)) == NULL)
        XMLException::selfThrow("StrategyXML - ERROR : at least one " + ONESTEPINTEGRATOR_TAG + " must be declared.");
      while (OSInode != NULL)
      {
        type = (char*)OSInode->name;
        if (type == MOREAU_TAG)
          oneStepIntegratorXMLVector.push_back(new MoreauXML(OSInode, DSAvailabilityMap));

        else if (type == LSODAR_TAG)
          oneStepIntegratorXMLVector.push_back(new LsodarXML(OSInode, DSAvailabilityMap));
        else
          XMLException::selfThrow("StrategyXML, Integrator loading : undefined OneStepIntegrator type : " + type);

        OSInode = SiconosDOMTreeTools::findFollowNode(OSInode);
      }
    }
    else
      XMLException::selfThrow("StrategyXML - ERROR : tag " + ONESTEPINTEGRATOR_DEFINITION_TAG + " not found.");
  }
  else cout << "StrategyXML - strategy XML constructor  : no integrators defined, use of LMGC90 tag." << endl;

  // === OneStepNSProblem data loading ===
  if ((node = SiconosDOMTreeTools::findNodeChild(strategyNode, ONESTEPNSPROBLEM_TAG)) != NULL)
  {
    xmlNode *NSnode = SiconosDOMTreeTools::findNodeChild(node);

    // !!! first node MUST be formulation for solving ie NSPB type!!!
    if (NSnode != NULL)
    {
      string type((char*)NSnode->name);
      if (type == LCP_TAG)
        oneStepNSProblemXML = new LCPXML(node, definedNumberInteractionVector);

      else if (type == QP_TAG)
        oneStepNSProblemXML = new QPXML(node, definedNumberInteractionVector);

      else if (type == FrictionContact2D_TAG || type == FrictionContact3D_TAG)
        oneStepNSProblemXML = new FrictionContactXML(node, definedNumberInteractionVector);

      else if (type == RELAY_TAG) //--Not implemented for the moment
        XMLException::selfThrow("StrategyXML constructor, following OneStepNSProblem is not yet implemented: " + type);

      else
        XMLException::selfThrow("StrategyXML constructor, undefined OneStepNSProblem type : " + type);
    }
  }
}


StrategyXML::~StrategyXML()
{
  if (timeDiscretisationXML != NULL)
    delete timeDiscretisationXML;

  if (oneStepNSProblemXML != NULL)
    delete oneStepNSProblemXML;

  if (oneStepIntegratorXMLVector.size() > 0)
  {
    for (unsigned int i = 0; i < oneStepIntegratorXMLVector.size(); i++)
      delete oneStepIntegratorXMLVector[i];
    oneStepIntegratorXMLVector.clear();
  }
}

// To be reviewed ...
void StrategyXML::saveStrategy2XML(xmlNode* node, Strategy* str)
{
  strategyNode = node;
  string type, tmp;
  xmlNode *integratorDefinitionNode;
  OneStepIntegratorXML* osixml;
  OneStepNSProblemXML* osnspbxml;
  TimeDiscretisationXML* tdxml;
  int i;

  if (strategyNode != NULL)
  {
    // === TimeDiscretisation node ===
    if (str->getTimeDiscretisationPtr()->getTimeDiscretisationXMLPtr() == NULL)
    {
      node = xmlNewChild(strategyNode, NULL, (xmlChar*)TIMEDISCRETISATION_TAG.c_str(), NULL);
      if (str->getTimeDiscretisationPtr()->isConstant())
        xmlNewProp(node, (xmlChar*)TD_ISCONSTANT.c_str(), (xmlChar*)"true");

      tdxml = new TimeDiscretisationXML();

      // linkage between the TimeDiscretisation and his TimeDiscretisationXML
      str->getTimeDiscretisationPtr()->setTimeDiscretisationXMLPtr(tdxml);

      // creation of the TimeDiscretisationXML
      tdxml->updateTimeDiscretisationXML(node, str->getTimeDiscretisationPtr());

      timeDiscretisationXML = tdxml;
    }

    if (SiconosDOMTreeTools::findNodeChild((const xmlNode*)strategyNode, LMGC90_STRATEGY_TAG) == NULL)
    {
      // === integrator_Definition node ===
      if (!hasOneStepIntegratorXML())
        integratorDefinitionNode = xmlNewChild(strategyNode, NULL, (xmlChar*)ONESTEPINTEGRATOR_DEFINITION_TAG.c_str(), NULL);

      // creation of the OneStepIntegratorXML objects
      for (i = 0; i < str->getOneStepIntegratorVectorSize(); i++)
      {
        if (str->getOneStepIntegrator(i)->getOneStepIntegratorXMLPtr() == NULL)
        {
          type = str->getOneStepIntegrator(i)->getType();
          if (type == MOREAU_TAG)
          {
            node = xmlNewChild(integratorDefinitionNode, NULL, (xmlChar*)MOREAU_TAG.c_str(), NULL);
            osixml = new MoreauXML();

            // linkage between the OneStepIntegrator and his OneStepIntegratorXML
            str->getOneStepIntegrator(i)->setOneStepIntegratorXMLPtr(osixml);

            // creation of the OneStepIntegratorXML
            static_cast<MoreauXML*>(osixml)->updateOneStepIntegratorXML(node, str->getOneStepIntegrator(i));

            oneStepIntegratorXMLVector.push_back(osixml);
          }
          else if (type == LSODAR_TAG)
          {
            node = xmlNewChild(integratorDefinitionNode, NULL, (xmlChar*)LSODAR_TAG.c_str(), NULL);
            osixml = new LsodarXML();

            // linkage between the OneStepIntegrator and his OneStepIntegratorXML
            str->getOneStepIntegrator(i)->setOneStepIntegratorXMLPtr(osixml);

            // creation of the OneStepIntegratorXML
            static_cast<LsodarXML*>(osixml)->updateOneStepIntegratorXML(node, str->getOneStepIntegrator(i));

            this->oneStepIntegratorXMLVector.push_back(osixml);
          }
          else
            XMLException::selfThrow("StrategyXML - saveStrategy2XML ERROR : undefined integrator type : " + type);
        }
      }
    }
    else
      XMLException::selfThrow("StrategyXML - saveStrategy2XML ERROR : LMGC90 strategy not yet implemented");

    // === OneStepNSProblemXML ===
    if (str->getOneStepNSProblemPtr() != NULL)
    {
      if (!hasOneStepNSProblemXML())
        node = xmlNewChild(strategyNode, NULL, (xmlChar*)ONESTEPNSPROBLEM_TAG.c_str(), NULL);

      if (str->getOneStepNSProblemPtr()->getOneStepNSProblemXML() == NULL)
      {
        type = str->getOneStepNSProblemPtr()->getType();
        if (type == LCP_TAG)
        {
          xmlNewChild(node, NULL, (xmlChar*)LCP_TAG.c_str(), NULL);
          osnspbxml = new LCPXML();

          // linkage between the OneStepNSProblem and his OneStepNSProblemXML
          str->getOneStepNSProblemPtr()->setOneStepNSProblemXML(osnspbxml);

          // creation of the OneStepNSProblemXML
          osnspbxml->updateOneStepNSProblemXML(node, str->getOneStepNSProblemPtr());

          oneStepNSProblemXML = osnspbxml;
        }
        else if (type == QP_TAG)
        {
          xmlNewChild(node, NULL, (xmlChar*)QP_TAG.c_str(), NULL);
          osnspbxml = new QPXML();

          // linkage between the OneStepNSProblem and his OneStepNSProblemXML
          str->getOneStepNSProblemPtr()->setOneStepNSProblemXML(osnspbxml);

          // creation of the OneStepNSProblemXML
          osnspbxml->updateOneStepNSProblemXML(node, str->getOneStepNSProblemPtr());

          oneStepNSProblemXML = osnspbxml;
        }
        else
          XMLException::selfThrow("StrategyXML - loadStrategy ERROR : undefined OneStepNSProblem type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
      }
    }
  }
  else XMLException::selfThrow("StrategyXML - loadStrategy ERROR : no strategyNode defined.");
  OUT("StrategyXML::loadStrategy\n");
}

