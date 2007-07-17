/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
#include "SimulationXML.h"

#include "OneStepIntegratorXML.h"
#include "OneStepNSProblemXML.h"
#include "TimeDiscretisationXML.h"

// includes to be deleted thanks to factories ?

#include "LsodarXML.h"
#include "MoreauXML.h"
#include "LCPXML.h"
#include "QPXML.h"
#include "FrictionContactXML.h"

using namespace std;

SimulationXML::SimulationXML():
  rootNode(NULL), oneStepNSProblemXML(NULL), timeDiscretisationXML(NULL)
{}

SimulationXML::SimulationXML(xmlNodePtr rootSimulationNode):
  rootNode(rootSimulationNode), oneStepNSProblemXML(NULL), timeDiscretisationXML(NULL)
{
  xmlNodePtr node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LMGC90_SIMULATION_TAG)) == NULL)
  {
    // === TimeDiscretisation data loading ===
    if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, TIMEDISCRETISATION_TAG)) != NULL)
      timeDiscretisationXML = new TimeDiscretisationXML(node);
    else
      XMLException::selfThrow("SimulationXML - simulation XML constructor  ERROR : tag " + TIMEDISCRETISATION_TAG + " not found.");

    // === OneStepIntegrator data loading ===
    if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, ONESTEPINTEGRATOR_DEFINITION_TAG)) != NULL)
    {
      xmlNodePtr OSINode = SiconosDOMTreeTools::findNodeChild(node);
      if (OSINode == NULL) // At least one OSI must be described in the xml file.
        XMLException::selfThrow("SimulationXML - ERROR : at least one " + ONESTEPINTEGRATOR_TAG + " must be declared.");

      string typeOSI; // OneStepIntegrator type
      while (OSINode != NULL)
      {
        typeOSI = (char*)OSINode->name;
        if (typeOSI == MOREAU_TAG)
          OSIXMLSet.insert(new MoreauXML(OSINode));

        else if (typeOSI == LSODAR_TAG)
          OSIXMLSet.insert(new LsodarXML(OSINode));

        else
          XMLException::selfThrow("SimulationXML, Integrator loading : undefined OneStepIntegrator type: " + typeOSI);

        // go to next node
        OSINode = SiconosDOMTreeTools::findFollowNode(OSINode);
      }
    }
    else
      XMLException::selfThrow("SimulationXML - ERROR : tag " + ONESTEPINTEGRATOR_DEFINITION_TAG + " not found.");
  }
  else cout << "SimulationXML - Constructor : the Simulation is not defined -> the LMGC90 tag is used." << endl;

  // === OneStepNSProblem data loading ===
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, ONESTEPNSPROBLEM_TAG)) != NULL)
  {
    xmlNodePtr NSnode = SiconosDOMTreeTools::findNodeChild(node);

    // !!! first node MUST be formulation for solving ie NSPB type!!!
    if (NSnode != NULL)
    {
      string type((char*)NSnode->name);
      if (type == LCP_TAG)
        oneStepNSProblemXML = new LCPXML(node);

      else if (type == QP_TAG)
        oneStepNSProblemXML = new QPXML(node);

      else if (type == FrictionContact2D_TAG || type == FrictionContact3D_TAG)
        oneStepNSProblemXML = new FrictionContactXML(node);

      else if (type == RELAY_TAG) //--Not implemented for the moment
        XMLException::selfThrow("SimulationXML constructor, following OneStepNSProblem is not yet implemented: " + type);

      else
        XMLException::selfThrow("SimulationXML constructor, undefined OneStepNSProblem type : " + type);
    }
  }
}


SimulationXML::~SimulationXML()
{
  if (timeDiscretisationXML != NULL)
    delete timeDiscretisationXML;

  if (oneStepNSProblemXML != NULL)
    delete oneStepNSProblemXML;

  // Delete OSIXML set ...
  SetOfOSIXMLIt it;
  for (it = OSIXMLSet.begin(); it != OSIXMLSet.end(); ++it)
    if ((*it) != NULL) delete(*it);
  OSIXMLSet.clear();
}

// To be reviewed ...
void SimulationXML::saveSimulation2XML(xmlNodePtr  node, Simulation* str)
{
  XMLException::selfThrow("SimulationXML saveSimulation2XML, not yet implemented");

  //   rootNode = node;
  //   string type, tmp;
  //   xmlNodePtr integratorDefinitionNode;
  //   OneStepIntegratorXML* osixml;
  //   OneStepNSProblemXML* osnspbxml;
  //   TimeDiscretisationXML* tdxml;
  //   int i;

  //   if( rootNode != NULL )
  //     {
  //       // === TimeDiscretisation node ===
  //       if( str->getTimeDiscretisationPtr()->getTimeDiscretisationXMLPtr() == NULL )
  //  {
  //    node = xmlNewChild( rootNode, NULL, (xmlChar*)TIMEDISCRETISATION_TAG.c_str(), NULL );
  //    if(str->getTimeDiscretisationPtr()->isConstant())
  //      xmlNewProp( node, (xmlChar*)TD_ISCONSTANT.c_str(), (xmlChar*)"true" );

  //    tdxml = new TimeDiscretisationXML();

  //    // linkage between the TimeDiscretisation and his TimeDiscretisationXML
  //    str->getTimeDiscretisationPtr()->setTimeDiscretisationXMLPtr( tdxml );

  //    // creation of the TimeDiscretisationXML
  //    tdxml->updateTimeDiscretisationXML( node, str->getTimeDiscretisationPtr() );

  //    timeDiscretisationXML = tdxml;
  //  }

  //       if( SiconosDOMTreeTools::findNodeChild((const xmlNodePtr )rootNode, LMGC90_SIMULATION_TAG) == NULL )
  //  {
  //    // === integrator_Definition node ===
  //    if( !hasOneStepIntegratorXML() )
  //      integratorDefinitionNode = xmlNewChild(rootNode, NULL, (xmlChar*)ONESTEPINTEGRATOR_DEFINITION_TAG.c_str(), NULL);

  //    // creation of the OneStepIntegratorXML objects
  //    for(i=0; i < str->getOneStepIntegratorVectorSize(); i++)
  //      {
  //        if( str->getOneStepIntegrator(i)->getOneStepIntegratorXMLPtr() == NULL )
  //    {
  //      type = str->getOneStepIntegrator(i)->getType();
  //      if (type == MOREAU_TAG)
  //        {
  //          node = xmlNewChild( integratorDefinitionNode, NULL, (xmlChar*)MOREAU_TAG.c_str(), NULL );
  //          osixml = new MoreauXML();

  //          // linkage between the OneStepIntegrator and his OneStepIntegratorXML
  //          str->getOneStepIntegrator(i)->setOneStepIntegratorXMLPtr( osixml );

  //          // creation of the OneStepIntegratorXML
  //          static_cast<MoreauXML*>(osixml)->updateOneStepIntegratorXML( node, str->getOneStepIntegrator(i) );

  //          oneStepIntegratorXMLVector.push_back( osixml );
  //        }
  //      else if (type == LSODAR_TAG)
  //        {
  //          node = xmlNewChild( integratorDefinitionNode, NULL, (xmlChar*)LSODAR_TAG.c_str(), NULL );
  //          osixml = new LsodarXML();

  //          // linkage between the OneStepIntegrator and his OneStepIntegratorXML
  //          str->getOneStepIntegrator(i)->setOneStepIntegratorXMLPtr( osixml );

  //          // creation of the OneStepIntegratorXML
  //          static_cast<LsodarXML*>(osixml)->updateOneStepIntegratorXML( node, str->getOneStepIntegrator(i) );

  //          this->oneStepIntegratorXMLVector.push_back( osixml );
  //        }
  //      else
  //        XMLException::selfThrow("SimulationXML - saveSimulation2XML ERROR : undefined integrator type : " + type );
  //    }
  //      }
  //  }
  //       else
  //  XMLException::selfThrow("SimulationXML - saveSimulation2XML ERROR : LMGC90 simulation not yet implemented");

  //       // === OneStepNSProblemXML ===
  //       if( str->getOneStepNSProblemPtr() != NULL )
  //  {
  //    if( !hasOneStepNSProblemXML() )
  //      node = xmlNewChild(rootNode, NULL, (xmlChar*)ONESTEPNSPROBLEM_TAG.c_str(), NULL );

  //    if( str->getOneStepNSProblemPtr()->getOneStepNSProblemXML() == NULL )
  //      {
  //        type = str->getOneStepNSProblemPtr()->getType();
  //        if (type == LCP_TAG)
  //    {
  //      xmlNewChild( node, NULL, (xmlChar*)LCP_TAG.c_str(), NULL );
  //      osnspbxml = new LCPXML();

  //      // linkage between the OneStepNSProblem and his OneStepNSProblemXML
  //      str->getOneStepNSProblemPtr()->setOneStepNSProblemXML( osnspbxml );

  //      // creation of the OneStepNSProblemXML
  //      osnspbxml->updateOneStepNSProblemXML( node, str->getOneStepNSProblemPtr() );

  //      oneStepNSProblemXML = osnspbxml;
  //    }
  //        else if (type == QP_TAG)
  //    {
  //      xmlNewChild( node, NULL, (xmlChar*)QP_TAG.c_str(), NULL );
  //      osnspbxml = new QPXML();

  //      // linkage between the OneStepNSProblem and his OneStepNSProblemXML
  //      str->getOneStepNSProblemPtr()->setOneStepNSProblemXML( osnspbxml );

  //      // creation of the OneStepNSProblemXML
  //      osnspbxml->updateOneStepNSProblemXML( node, str->getOneStepNSProblemPtr() );

  //      oneStepNSProblemXML = osnspbxml;
  //    }
  //        else
  //    XMLException::selfThrow("SimulationXML - loadSimulation ERROR : undefined OneStepNSProblem type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
  //      }
  //  }
  //     }
  //   else XMLException::selfThrow("SimulationXML - loadSimulation ERROR : no rootNode defined.");
}

