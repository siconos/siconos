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

#include "SiconosPointers.hpp"
#include "SimulationXML.hpp"

#include "OneStepIntegratorXML.hpp"
#include "OneStepNSProblemXML.hpp"
#include "TimeDiscretisationXML.hpp"

// includes to be deleted thanks to factories ?

#include "LsodarXML.hpp"
#include "MoreauXML.hpp"
#include "QPXML.hpp"
#include "FrictionContactXML.hpp"

using namespace std;

SimulationXML::SimulationXML(xmlNodePtr rootSimulationNode): rootNode(rootSimulationNode)
{
  xmlNodePtr node;

  if (!(node = SiconosDOMTreeTools::findNodeChild(rootNode, LMGC90_SIMULATION_TAG)))
  {
    // === TimeDiscretisation data loading ===
    if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, TIMEDISCRETISATION_TAG)))
      _timeDiscretisationXML.reset(new TimeDiscretisationXML(node));
    else
      XMLException::selfThrow("SimulationXML - simulation XML constructor  ERROR : tag " + TIMEDISCRETISATION_TAG + " not found.");

    // === OneStepIntegrator data loading ===
    if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, ONESTEPINTEGRATOR_DEFINITION_TAG)))
    {
      xmlNodePtr OSINode = SiconosDOMTreeTools::findNodeChild(node);
      if (!OSINode)  // At least one OSI must be described in the xml file.
        XMLException::selfThrow("SimulationXML - ERROR : at least one " + ONESTEPINTEGRATOR_TAG + " must be declared.");

      string typeOSI; // OneStepIntegrator type
      while (OSINode)
      {
        typeOSI = (char*)OSINode->name;
        if (typeOSI == MOREAU_TAG)
          OSIXMLSet.insert(SP::MoreauXML(new MoreauXML(OSINode)));

        else if (typeOSI == LSODAR_TAG)
          OSIXMLSet.insert(SP::LsodarXML(new LsodarXML(OSINode)));

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
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, ONESTEPNSPROBLEM_TAG)))
  {
    xmlNodePtr OSNSPBNode = SiconosDOMTreeTools::findNodeChild(node);

    string typeOSNS; // OneStepNSPb type
    while (OSNSPBNode)
    {
      typeOSNS = (char*)OSNSPBNode->name;
      if (typeOSNS == LCP_TAG) // No need to specific LCPXML
        OSNSPBXMLSet.insert(SP::OneStepNSProblemXML(new OneStepNSProblemXML(OSNSPBNode)));

      else if (typeOSNS == QP_TAG)
        OSNSPBXMLSet.insert(SP::QPXML(new QPXML(OSNSPBNode)));

      else if (typeOSNS == "FrictionContact" || typeOSNS == "PrimalFrictionContact")
        OSNSPBXMLSet.insert(SP::FrictionContactXML(new FrictionContactXML(OSNSPBNode)));

      else // if (typeOSNS == RELAY_TAG) //--Not implemented for the moment
        XMLException::selfThrow("SimulationXML constructor, undefined (or not yet implemented) OneStepNSProblem of type: " + typeOSNS);

      // go to next node
      OSNSPBNode = SiconosDOMTreeTools::findFollowNode(OSNSPBNode);
    }
  }
  // Else nothing: it is possible to define a simulation without any OneStepNSProblem
}

SimulationXML::~SimulationXML()
{
  OSIXMLSet.clear();
  OSNSPBXMLSet.clear();
}

// To be reviewed ...
void SimulationXML::saveSimulation2XML(xmlNodePtr  node, SP::Simulation str)
{
  XMLException::selfThrow("SimulationXML saveSimulation2XML, not yet implemented");

  //   rootNode = node;
  //   string type, tmp;
  //   xmlNodePtr integratorDefinitionNode;
  //   OneStepIntegratorXML* osixml;
  //   SP::OneStepNSProblemXML  osnspbxml;
  //   TimeDiscretisationXML* tdxml;
  //   int i;

  //   if( rootNode != NULL )
  //     {
  //       // === TimeDiscretisation node ===
  //       if( str->timeDiscretisation()->timeDiscretisationXML() == NULL )
  //  {
  //    node = xmlNewChild( rootNode, NULL, (xmlChar*)TIMEDISCRETISATION_TAG.c_str(), NULL );
  //    if(str->timeDiscretisation()->isConstant())
  //      xmlNewProp( node, (xmlChar*)TD_ISCONSTANT.c_str(), (xmlChar*)"true" );

  //    tdxml = new TimeDiscretisationXML();

  //    // linkage between the TimeDiscretisation and his TimeDiscretisationXML
  //    str->timeDiscretisation()->setTimeDiscretisationXMLPtr( tdxml );

  //    // creation of the TimeDiscretisationXML
  //    tdxml->updateTimeDiscretisationXML( node, str->timeDiscretisation() );

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
  //        if( str->getOneStepIntegrator(i)->oneStepIntegratorXML() == NULL )
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
  //       if( str->oneStepNSProblem() != NULL )
  //  {
  //    if( !hasOneStepNSProblemXML() )
  //      node = xmlNewChild(rootNode, NULL, (xmlChar*)ONESTEPNSPROBLEM_TAG.c_str(), NULL );

  //    if( str->oneStepNSProblem()->getOneStepNSProblemXML() == NULL )
  //      {
  //        type = str->oneStepNSProblem()->getType();
  //        if (type == LCP_TAG)
  //    {
  //      xmlNewChild( node, NULL, (xmlChar*)LCP_TAG.c_str(), NULL );
  //      osnspbxml = new LCPXML();

  //      // linkage between the OneStepNSProblem and his OneStepNSProblemXML
  //      str->oneStepNSProblem()->setOneStepNSProblemXML( osnspbxml );

  //      // creation of the OneStepNSProblemXML
  //      osnspbxml->updateOneStepNSProblemXML( node, str->oneStepNSProblem() );

  //      oneStepNSProblemXML = osnspbxml;
  //    }
  //        else if (type == QP_TAG)
  //    {
  //      xmlNewChild( node, NULL, (xmlChar*)QP_TAG.c_str(), NULL );
  //      osnspbxml = new QPXML();

  //      // linkage between the OneStepNSProblem and his OneStepNSProblemXML
  //      str->oneStepNSProblem()->setOneStepNSProblemXML( osnspbxml );

  //      // creation of the OneStepNSProblemXML
  //      osnspbxml->updateOneStepNSProblemXML( node, str->oneStepNSProblem() );

  //      oneStepNSProblemXML = osnspbxml;
  //    }
  //        else
  //    XMLException::selfThrow("SimulationXML - loadSimulation ERROR : undefined OneStepNSProblem type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
  //      }
  //  }
  //     }
  //   else XMLException::selfThrow("SimulationXML - loadSimulation ERROR : no rootNode defined.");
}

