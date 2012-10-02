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
#include "NonSmoothDynamicalSystemXML.hpp"

//---  Following includes to be suppressed thanks to factory (?) ---

// DS
#include "LagrangianDSXML.hpp"
#include "InteractionXML.hpp"
#include "LagrangianLinearTIDSXML.hpp"
#include "FirstOrderLinearDSXML.hpp"

using namespace std;


NonSmoothDynamicalSystemXML::NonSmoothDynamicalSystemXML(): rootNode(NULL)
{}

NonSmoothDynamicalSystemXML::NonSmoothDynamicalSystemXML(xmlNodePtr  rootNSDSNode): rootNode(rootNSDSNode)
{
  if (rootNode)
  {
    xmlNodePtr node;

    if (!(node = SiconosDOMTreeTools::findNodeChild(rootNode, LMGC90_NSDS_TAG)))
    {
      if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, DYNAMICAL_SYSTEM_DEFINITION_TAG)))
        loadDynamicalSystemXML(node);
      else
        XMLException::selfThrow("NonSmoothDynamicalSystemXML - loadNonSmoothDynamicalSystem error : tag " + DYNAMICAL_SYSTEM_DEFINITION_TAG + " not found.");

      // === Interactions ===
      if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, INTERACTION_DEFINITION_TAG)))
      {
        xmlNodePtr interNode = SiconosDOMTreeTools::findNodeChild(node, INTERACTION_TAG);
        // We look for all Interaction tag, and for each of them add an InteractionXML pointer in the interactionXMLSet
        while (interNode) // scan all the "Interaction" tags, and for each of them insert an InteractionXML object in the set
        {
          interactionsXMLSet.insert(SP::InteractionXML(new InteractionXML(interNode)));
          interNode = SiconosDOMTreeTools::findFollowNode(interNode, INTERACTION_TAG);
        }
      }
    }
    else cout << "NonSmoothDynamicalSystemXML -Constructor: the Non Smooth Dynamical System is not defined -> the LMGC90 tag is used." << endl;
  }
}

NonSmoothDynamicalSystemXML::~NonSmoothDynamicalSystemXML()
{
  // Delete DSXML set ...
  DSXMLSet.clear();

  // Delete InteractionXML set ...
  interactionsXMLSet.clear();
}

void NonSmoothDynamicalSystemXML::loadDynamicalSystemXML(xmlNodePtr  rootDSNode)
{
  xmlNodePtr node;

  string type; //Type of DS
  bool isbvp = isBVP();

  // rootDSNode = "DS_Definition". We look for its children node (DynamicalSystem and derived classes) and for
  // each of them add a new DynamicalSystemXML in the set of DSXML.
  node = SiconosDOMTreeTools::findNodeChild((xmlNodePtr)rootDSNode);
  if (!node)  // At least one DynamicalSystem must be described in the xml file.
    XMLException::selfThrow("NonSmoothDynamicalSystemXML - loadDynamicalSystemXML error : at least one " + DYNAMICAL_SYSTEM_TAG + " must be declared.");

  while (node)
  {
    type = (char*)node->name; // get the type of DS
    if (type == LAGRANGIAN_NON_LINEARDS_TAG)
      DSXMLSet.insert(SP::LagrangianDSXML(new LagrangianDSXML(node, isbvp)));
    else if (type == LAGRANGIAN_TIDS_TAG)
      DSXMLSet.insert(SP::LagrangianLinearTIDSXML(new LagrangianLinearTIDSXML(node, isbvp)));
    else if (type == LINEAR_DS_TAG || type == LINEAR_TIDS_TAG)
      DSXMLSet.insert(SP::FirstOrderLinearDSXML(new FirstOrderLinearDSXML(node, isbvp)));
    else if (type == NON_LINEAR_DS_TAG)
      DSXMLSet.insert(SP::FirstOrderNonLinearDSXML(new FirstOrderNonLinearDSXML(node, isbvp)));
    else
      XMLException::selfThrow("NonSmoothDynamicalSystemXML - loadDynamicalSystemXML error : undefined DS type: " + type);
    // go to next node ...
    node = SiconosDOMTreeTools::findFollowNode(node);
  }
}

// WARNING : FOLLOWING FUNCTIONS ARE OBSOLETE, USELESS OR AT LEAST TO BE REVIEWED WHEN DSIO AND EQUALITY CONSTRAINTS
// WILL BE WELL IMPLEMENTED IN MODELING PACKAGE

void NonSmoothDynamicalSystemXML::loadNonSmoothDynamicalSystem(SP::NonSmoothDynamicalSystem nsds)
{
  XMLException::selfThrow("NonSmoothDynamicalSystemXML - loadDynamicalSystemXML not implemented.");
  //   string type;
  //   string tmp;
  //   xmlNodePtr  node, ecDsioNode;
  //   xmlNodePtr  dsDefinitionNode;
  //   xmlNodePtr  interactionDefinitionNode, ecDefinitionNode;
  //   DynamicalSystemXML* dsxml;
  //   InteractionXML* interactionXML;
  //   EqualityConstraintXML *ecXML;
  //   int number;
  //   unsigned int i;
  //   char num[32];
  //   map<int, DynamicalSystemXML*>::iterator it;
  //   map<int, InteractionXML*>::iterator itinter;
  //   map<int, EqualityConstraintXML*>::iterator itec;

  //   if( rootNode )
  //     {
  //       setBVP( nsds->isBVP() );

  //       // at first, we check whether we the tag is LMGC90 tag
  //       if(! SiconosDOMTreeTools::findNodeChild((const xmlNodePtr )rootNode, LMGC90_NSDS_TAG)  )
  //  {
  //    // creation of the DS_Definition node if necessary
  //    dsDefinitionNode = SiconosDOMTreeTools::findNodeChild((const xmlNodePtr )rootNode, DYNAMICAL_SYSTEM_DEFINITION_TAG);
  //    if( !dsDefinitionNode )
  //      dsDefinitionNode = xmlNewChild(rootNode, NULL, (xmlChar*)DYNAMICAL_SYSTEM_DEFINITION_TAG.c_str(), NULL);

  //    /*
  //     * now, creation of the DynamicalSystemXML objects
  //     */
  //    for(i=0; i<nsds->getNumberOfDS(); i++)
  //      {
  //        if( !nsds->dynamicalSystem(i)->dynamicalSystemXML() )
  //    {
  //      type = nsds->dynamicalSystem(i)->getType();
  //      number = nsds->dynamicalSystem(i)->getNumber();
  //      sprintf(num, "%d", number);
  //      definedDSNumbers.push_back( number );


  //      // verifies if this Dynamical System has a number which not used
  //      it = DSXMLMap.find(number);
  //      if( it == DSXMLMap.end() )
  //        {
  //          //node = xmlNewChild( dsDefinitionNode, NULL, (xmlChar*)NSDS_DS.c_str(), NULL );
  //          if (type == Type::LagrangianDS)
  //      {
  //        node = xmlNewChild( dsDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_NON_LINEARDS_TAG.c_str(), NULL );
  //        xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //        dsxml = new LagrangianDSXML();

  //        // linkage between the DynamicalSystem and his DynamicalSystemXML
  //        nsds->dynamicalSystem(i)->setDynamicalSystemXMLPtr( dsxml );

  //        // creation of the DynamicalSystemXML
  //        static_cast<LagrangianDSXML*>(dsxml)->updateDynamicalSystemXML( node, nsds->dynamicalSystem(i) );

  //        DSXMLMap[number] = dsxml;
  //      }
  //          else if (type == LTIDS)
  //      {
  //        node = xmlNewChild( dsDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_TIDS_TAG.c_str(), NULL );
  //        xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //        dsxml = new LagrangianLinearTIDSXML();

  //        // linkage between the DynamicalSystem and his DynamicalSystemXML
  //        nsds->dynamicalSystem(i)->setDynamicalSystemXMLPtr( dsxml );

  //        // creation of the DynamicalSystemXML
  //        static_cast<LagrangianLinearTIDSXML*>(dsxml)->updateDynamicalSystemXML( node, nsds->dynamicalSystem(i) );

  //        DSXMLMap[number] = dsxml;
  //      }
  //          else if (type == LDS)
  //      {
  //        node = xmlNewChild( dsDefinitionNode, NULL, (xmlChar*)LINEAR_DS_TAG.c_str(), NULL );
  //        xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //        dsxml = new LinearDSXML();

  //        // linkage between the DynamicalSystem and his DynamicalSystemXML
  //        nsds->dynamicalSystem(i)->setDynamicalSystemXMLPtr( dsxml );

  //        // creation of the DynamicalSystemXML
  //        static_cast<LinearDSXML*>(dsxml)->updateDynamicalSystemXML( node, nsds->dynamicalSystem(i) );

  //        DSXMLMap[number] = dsxml;

  //      }
  //          else if (type == NLDS)
  //      {
  //        node = xmlNewChild( dsDefinitionNode, NULL, (xmlChar*)NON_LINEAR_DS_TAG.c_str(), NULL );
  //        xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //        dsxml = new DynamicalSystemXML();

  //        // linkage between the DynamicalSystem and his DynamicalSystemXML
  //        nsds->dynamicalSystem(i)->setDynamicalSystemXMLPtr( dsxml );

  //        // creation of the DynamicalSystemXML
  //        dsxml->updateDynamicalSystemXML( node, nsds->dynamicalSystem(i) );

  //        DSXMLMap[number] = dsxml;
  //      }
  //          else
  //      {
  //        XMLException::selfThrow("NonSmoothDynamicalSystemXML - loadNonSmoothDynamicalSystem error : undefined DS type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
  //      }
  //        }
  //      else
  //        {
  //          tmp = num;
  //          XMLException::selfThrow("NonSmoothDynamicalSystemXML - loadNonSmoothDynamicalSystem | Error : the Dynamical System number : " + tmp + " already exists!");
  //        }
  //    }
  //        else
  //    cout<<"## /!\\ the DynamicalSystem : "<<nsds->dynamicalSystem(i)->getType()<<" number "<<nsds->dynamicalSystem(i)->getNumber()<<
  //      ", has already an XML object."<<endl;
  //      }
  //  }
  //       else
  //  {
  //    // the LMGC90 tag for DS definition is in the XML file
  //    // \todo !!!!!!  => specific treatments todo
  //  }


  //       // creation of the EqualityConstraint_Defintion if necessary
  //       if( nsds->getEqualityConstraints().size() > 0 )
  //  {
  //    ecDefinitionNode = SiconosDOMTreeTools::findNodeChild((const xmlNodePtr )rootNode, EQUALITYCONSTRAINT_DEFINITION_TAG);
  //    if( !ecDefinitionNode )
  //      ecDefinitionNode = xmlNewChild(rootNode, NULL, (xmlChar*)EQUALITYCONSTRAINT_DEFINITION_TAG.c_str(), NULL);

  //    for(i=0; i<nsds->getEqualityConstraints().size(); i++)
  //      {
  //        if(! nsds->equalityConstraint(i)->getEqualityConstraintXML() )
  //    {
  //      number = nsds->equalityConstraint(i)->getNumber();
  //      sprintf(num, "%d", number);
  //      definedEqualityConstraintNumbers.push_back( number );

  //      //verifies if the EqualityConstraint has been defined before
  //      itec = equalityConstraintXMLMap.find(number);
  //      if (itec == equalityConstraintXMLMap.end())
  //        {
  //          if( nsds->equalityConstraint(i)->getType() == LINEAREC )
  //      {
  //        node = xmlNewChild( ecDefinitionNode, NULL, (xmlChar*)LINEAR_EC_TAG.c_str(), NULL );
  //        xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //        ecXML = new LinearECXML();

  //        // linkage between the DynamicalSystem and his DynamicalSystemXML
  //        nsds->equalityConstraint(i)->setEqualityConstraintXML( ecXML );

  //        // creation of the DynamicalSystemXML
  //        static_cast<LinearECXML*>(ecXML)->updateEqualityConstraintXML( node, nsds->equalityConstraint(i) );

  //        equalityConstraintXMLMap[number] = ecXML;
  //      }
  //          else if( nsds->equalityConstraint(i)->getType() == LINEARTIEC )
  //      {
  //        node = xmlNewChild( ecDefinitionNode, NULL, (xmlChar*)LINEAR_TIME_INVARIANT_EC_TAG.c_str(), NULL );
  //        xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //        ecXML = new LinearTIECXML();

  //        // linkage between the DynamicalSystem and his DynamicalSystemXML
  //        nsds->equalityConstraint(i)->setEqualityConstraintXML( ecXML );

  //        // creation of the DynamicalSystemXML
  //        static_cast<LinearTIECXML*>(ecXML)->updateEqualityConstraintXML( node, nsds->equalityConstraint(i) );

  //        equalityConstraintXMLMap[number] = ecXML;
  //      }
  //          else if( nsds->equalityConstraint(i)->getType() == LAGRANGIANEC )
  //      {
  //        node = xmlNewChild( ecDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_EC_TAG.c_str(), NULL );
  //        xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //        ecXML = new LagrangianECXML();

  //        // linkage between the DynamicalSystem and his DynamicalSystemXML
  //        nsds->equalityConstraint(i)->setEqualityConstraintXML( ecXML );

  //        // creation of the DynamicalSystemXML
  //        static_cast<LagrangianECXML*>(ecXML)->updateEqualityConstraintXML( node, nsds->equalityConstraint(i) );

  //        equalityConstraintXMLMap[number] = ecXML;
  //      }
  //          else if( nsds->equalityConstraint(i)->getType() == LAGRANGIANLINEAREC )
  //      {
  //        node = xmlNewChild( ecDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_LINEAR_EC_TAG.c_str(), NULL );
  //        xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //        ecXML = new LagrangianECXML();

  //        // linkage between the DynamicalSystem and his DynamicalSystemXML
  //        nsds->equalityConstraint(i)->setEqualityConstraintXML( ecXML );

  //        // creation of the DynamicalSystemXML
  //        static_cast<LagrangianLinearECXML*>(ecXML)->updateEqualityConstraintXML( node, nsds->equalityConstraint(i) );

  //        equalityConstraintXMLMap[number] = ecXML;
  //      }
  //          else if( nsds->equalityConstraint(i)->getType() == NLINEAREC )
  //      {
  //        node = xmlNewChild( ecDefinitionNode, NULL, (xmlChar*)NON_LINEAR_EC_TAG.c_str(), NULL );
  //        xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //        ecXML = new EqualityConstraintXML();

  //        // linkage between the DynamicalSystem and his DynamicalSystemXML
  //        nsds->equalityConstraint(i)->setEqualityConstraintXML( ecXML );

  //        // creation of the DynamicalSystemXML
  //        ecXML->updateEqualityConstraintXML( node, nsds->equalityConstraint(i) );

  //        equalityConstraintXMLMap[number] = ecXML;
  //      }
  //          else XMLException::selfThrow("NonSmoothDynamicalSystemXML - loadNonSmoothDynamicalSystem | Error : the EqualityConstraint type : " + nsds->equalityConstraint(i)->getType() + " doesn't exist!");

  //          /*  end of the save : saving the DynamicalSystem linked to this DSInputOutput */
  //          ecDsioNode = xmlNewChild( node, NULL, (xmlChar*)DSIO_CONCERNED.c_str(), NULL );
  //          for( unsigned int j=0; j<nsds->equalityConstraint(i)->getDSInputOutputs().size(); j++)
  //      {
  //        node = xmlNewChild( ecDsioNode, NULL, (xmlChar*)DSINPUTOUTPUT_TAG.c_str(), NULL );
  //        number = nsds->equalityConstraint(i)->getDSInputOutput(j)->getNumber();
  //        sprintf(num, "%d", number);
  //        xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //      }
  //        }
  //      else
  //        {
  //          tmp = num;
  //          XMLException::selfThrow("NonSmoothDynamicalSystemXML - loadNonSmoothDynamicalSystem | Error : the EqualityConstraint number : " + tmp + " already exists!");
  //        }
  //    }
  //      }
  //  }


  //       // creation of the Interaction_Defintion if necessary
  //       if( nsds->getInteractionVectorSize() > 0 )
  //  {
  //    interactionDefinitionNode = SiconosDOMTreeTools::findNodeChild((const xmlNodePtr )rootNode, INTERACTION_DEFINITION_TAG);
  //    if( !interactionDefinitionNode )
  //      interactionDefinitionNode = xmlNewChild(rootNode, NULL, (xmlChar*)INTERACTION_DEFINITION_TAG.c_str(), NULL);

  //    for(i=0; int(i)<nsds->getInteractionVectorSize(); i++)
  //      {
  //        if( !nsds->interaction(i)->interactionXML()  )
  //    {
  //      number = nsds->interaction(i)->getNumber();
  //      sprintf(num, "%d", number);
  //      definedInteractionNumbers.push_back( number );

  //      // verifies if this Dynamical System has a number which not used
  //      itinter = interactionXMLMap.find(number);
  //      if (itinter == interactionXMLMap.end())
  //        {
  //          node = xmlNewChild( interactionDefinitionNode, NULL, (xmlChar*)INTERACTION_TAG.c_str(), NULL );
  //          xmlNewProp( node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num );
  //          interactionXML = new InteractionXML();

  //          // linkage between the DynamicalSystem and his DynamicalSystemXML
  //          nsds->interaction(i)->setInteractionXMLPtr( interactionXML );

  //          // creation of the DynamicalSystemXML
  //          interactionXML->updateInteractionXML( node, nsds->interaction(i) );

  //          interactionXMLMap[number] = interactionXML;
  //        }
  //      else
  //        {
  //          tmp = num;
  //          XMLException::selfThrow("NonSmoothDynamicalSystemXML - loadNonSmoothDynamicalSystem | Error : the Interaction number : " + tmp + " already exists!");
  //        }
  //    }
  //      }
  //  }
  //     }
  //   else XMLException::selfThrow("NonSmoothDynamicalSystemXML - loadNonSmoothDynamicalSystem( SP::NonSmoothDynamicalSystem nsds ) Error : no rootNode is defined.");
}

void NonSmoothDynamicalSystemXML::updateNonSmoothDynamicalSystemXML(xmlNodePtr  node, SP::NonSmoothDynamicalSystem nsds)
{
  rootNode = node;
  loadNonSmoothDynamicalSystem(nsds);
}

