/* Siconos version 1.0, Copyright INRIA 2005.
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
#include "InteractionXML.h"

//---  Following includes to be suppressed thanks to factory (?) ---
// Relations
#include "LinearTIRXML.h"
#include "LagrangianLinearRXML.h"
#include "LagrangianRXML.h"
// Nslaw
#include "ComplementarityConditionNSLXML.h"
#include "RelayNSLXML.h"
#include "NewtonImpactLawNSLXML.h"
#include "NewtonImpactFrictionNSLXML.h"

#include "NonSmoothLaw.h"
using namespace std;


InteractionXML::InteractionXML():
  rootInteractionXMLNode(NULL), idNode(NULL), nInterNode(NULL),
  yNode(NULL),   lambdaNode(NULL),
  dsConcernedNode(NULL), dsListNode(NULL),  relationXML(NULL), nSLawXML(NULL),
  isRelationXMLAllocatedIn(false), isNsLawXMLAllocatedIn(false)
{}

InteractionXML::InteractionXML(xmlNode * interactionNode, vector<int> definedNumbers):
  rootInteractionXMLNode(interactionNode), idNode(NULL), nInterNode(NULL),
  yNode(NULL),   lambdaNode(NULL),
  dsConcernedNode(NULL), dsListNode(NULL), relationXML(NULL), nSLawXML(NULL),
  isRelationXMLAllocatedIn(false), isNsLawXMLAllocatedIn(false)
{
  xmlNode *node, *node2;
  string type;

  // id (optional)
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, ID_ATTRIBUTE)) != NULL)
    idNode = node;
  // nInter (required)
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_NINTER)) != NULL)
    nInterNode = node;
  else
    XMLException::selfThrow("InteractionXML - loadInteractionProperties error : tag " + INTERACTION_NINTER + " not found.");

  // y (optional)
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_Y)) != NULL)
    yNode = node;
  // lambda (optional)
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_LAMBDA)) != NULL)
    lambdaNode = node;
  // dsConcerned (required)
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_DS_CONCERNED)) != NULL)
  {
    dsConcernedNode = node;
    // Check if all ds are concerned or not
    if (! hasAll())
    {
      // Get the DSList node
      if (SiconosDOMTreeTools::findNodeChild(dsConcernedNode, INDEX_LIST) != NULL)
        dsListNode = SiconosDOMTreeTools::findNodeChild(dsConcernedNode, INDEX_LIST);
      else
        XMLException::selfThrow("tag DSlist not found.");
    }
  }
  else
    XMLException::selfThrow("InteractionXML - xml constructor, tag " + INTERACTION_DS_CONCERNED + " not found.");

  // Relation and Non Smooth Law (required)
  if ((node = SiconosDOMTreeTools::findNodeChild((const xmlNode*)interactionNode, INTERACTION_CONTENT_TAG)) != NULL)
  {
    // the first child is the Relation
    if ((node2 = SiconosDOMTreeTools::findNodeChild(node)) != NULL)
    {
      // get Relation type
      type = (char*)node2->name;
      if (type == RELATION_TAG)
        relationXML = new RelationXML(node2);
      else if (type == LINEAR_TIME_INVARIANT_RELATION_TAG)
        relationXML = new LinearTIRXML(node2);
      else if (type == LAGRANGIAN_RELATION_TAG)
        relationXML = new LagrangianRXML(node2);
      else if (type == LAGRANGIAN_LINEAR_RELATION_TAG)
        relationXML = new LagrangianLinearRXML(node2);
      else
        XMLException::selfThrow("InteractionXML : undefined Relation type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?)");
      isRelationXMLAllocatedIn = true;
    }
    else
      XMLException::selfThrow("InteractionXML - loadInteractionProperties error : Relation not found.");

    // the second child is the NonSmoothLaw
    node = SiconosDOMTreeTools::findFollowNode(node2);
    if (node != NULL)
    {
      // Non smooth law type
      type = (char*)node->name;
      if (type == COMPLEMENTARITY_CONDITION_NSLAW_TAG)
        nSLawXML = new ComplementarityConditionNSLXML(node);
      else if (type == RELAY_NSLAW_TAG)
        nSLawXML = new RelayNSLXML(node);
      else if (type == NEWTON_IMPACT_LAW_NSLAW_TAG)
        nSLawXML = new NewtonImpactLawNSLXML(node);
      else if (type == NEWTON_IMPACT_FRICTION_NSLAW_TAG)
        nSLawXML = new NewtonImpactFrictionNSLXML(node);
      else
        XMLException::selfThrow("InteractionXML : undefined NonSmoothLaw type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?)");
      isNsLawXMLAllocatedIn = true;
    }
    else
      XMLException::selfThrow("InteractionXML - loadInteractionProperties error : NonSmoothLaw not found.");
  }
}

InteractionXML::~InteractionXML()
{
  if (isRelationXMLAllocatedIn)
  {
    delete relationXML;
    relationXML = NULL;
  }
  if (isNsLawXMLAllocatedIn)
  {
    delete nSLawXML;
    nSLawXML = NULL;
  }

}

void InteractionXML::updateInteractionXML(xmlNode* node, Interaction* inter)
{
  IN("InteractionXML::updateInteractionXML\n");
  rootInteractionXMLNode = node;
  loadInteraction(inter);
  OUT("InteractionXML::updateInteractionXML\n");
}

// warning: this function has not been reviewed for multiple DS loading !!
void InteractionXML::loadInteraction(Interaction* inter)
{
  IN("InteractionXML::loadInteraction( Interaction* )\n");
  string type;
  string tmp;
  xmlNode* node;
  RelationXML* newRelationXml;
  NonSmoothLawXML* nslxml;

  if (this->rootInteractionXMLNode != NULL)
  {
    /*
     * Creation of the RelationXML object
     */
    xmlNode *InteractionContentNode;
    if (inter->getRelationPtr() != NULL)
    {
      type = inter->getRelationPtr()->getType();
      InteractionContentNode = xmlNewChild(this->rootInteractionXMLNode, NULL, (xmlChar*)INTERACTION_CONTENT_TAG.c_str(), NULL);

      if (type == RELATION)
      {
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)RELATION_TAG.c_str(), NULL);
        relationXML = new RelationXML(node);
        // linkage between the Relation and his RelationXML
        inter->getRelationPtr()->setRelationXML(relationXML);
      }
      else if (type == LAGRANGIANLINEARRELATION)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_LL.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)LAGRANGIAN_LINEAR_RELATION_TAG.c_str(), NULL);
        newRelationXml = new LagrangianLinearRXML();

        // linkage between the Relation and his RelationXML
        inter->getRelationPtr()->setRelationXML(newRelationXml);

        // creation of the RelationXML
        static_cast<LagrangianLinearRXML*>(newRelationXml)->updateRelationXML(node, inter->getRelationPtr());

        relationXML = newRelationXml;
      }
      else if (type == LAGRANGIANRELATION)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_LNL.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)LAGRANGIAN_RELATION_TAG.c_str(), NULL);
        newRelationXml = new LagrangianRXML();

        // linkage between the Relation and his RelationXML
        inter->getRelationPtr()->setRelationXML(newRelationXml);

        // creation of the RelationXML
        static_cast<LagrangianRXML*>(newRelationXml)->updateRelationXML(node, inter->getRelationPtr());

        relationXML = newRelationXml;
      }
      else if (type == LINEARTIRELATION)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_LTI.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)LINEAR_TIME_INVARIANT_RELATION_TAG.c_str(), NULL);
        newRelationXml = new LinearTIRXML();

        // linkage between the Relation and his RelationXML
        inter->getRelationPtr()->setRelationXML(newRelationXml);

        // creation of the RelationXML
        static_cast<LinearTIRXML*>(newRelationXml)->updateRelationXML(node, inter->getRelationPtr());

        relationXML = newRelationXml;
      }
      else
      {
        XMLException::selfThrow("InteractionXML - loadInteraction error : undefined Relation type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
      }
    }
    else RuntimeException::selfThrow("InteractionXML::loadInteraction - There's no Relation in this Interaction, the XML platform can't be built");


    /*
     * now, creation of the NonSmoothLawXML object
     */
    if (inter->getNonSmoothLawPtr() != NULL)
    {
      type = inter->getNonSmoothLawPtr()->getType();
      //node = xmlNewChild( this->rootInteractionXMLNode, NULL, (xmlChar*)INTERACTION_NS_LAW.c_str(), NULL );
      if (type == COMPLEMENTARITYCONDITIONNSLAW)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_COMPLEMENTARITYCONDITIONNSLAW.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)COMPLEMENTARITY_CONDITION_NSLAW_TAG.c_str(), NULL);
        nslxml = new ComplementarityConditionNSLXML();

        // linkage between the Relation and his RelationXML
        inter->getNonSmoothLawPtr()->setNonSmoothLawXML(nslxml);

        // creation of the RelationXML
        static_cast<ComplementarityConditionNSLXML*>(nslxml)->updateNonSmoothLawXML(node, inter->getNonSmoothLawPtr());

        this->nSLawXML = nslxml;
      }
      else if (type == RELAYNSLAW)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_RELAYNSLAW.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)RELAY_NSLAW_TAG.c_str(), NULL);
        nslxml = new RelayNSLXML();

        // linkage between the Relation and his RelationXML
        inter->getNonSmoothLawPtr()->setNonSmoothLawXML(nslxml);

        // creation of the RelationXML
        static_cast<RelayNSLXML*>(nslxml)->updateNonSmoothLawXML(node, inter->getNonSmoothLawPtr());

        this->nSLawXML = nslxml;
      }
      else if (type == NEWTONIMPACTNSLAW)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_NEWTONIMPACTLAWNSLAW.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)NEWTON_IMPACT_LAW_NSLAW_TAG.c_str(), NULL);
        nslxml = new NewtonImpactLawNSLXML();

        // linkage between the Relation and his RelationXML
        inter->getNonSmoothLawPtr()->setNonSmoothLawXML(nslxml);

        // creation of the RelationXML
        static_cast<NewtonImpactLawNSLXML*>(nslxml)->updateNonSmoothLawXML(node, inter->getNonSmoothLawPtr());

        this->nSLawXML = nslxml;
      }
      else if (type == NEWTONIMPACTFRICTIONNSLAW)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_NEWTONIMPACTLAWNSLAW.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)NEWTON_IMPACT_FRICTION_NSLAW_TAG.c_str(), NULL);
        nslxml = new NewtonImpactFrictionNSLXML();

        // linkage between the Relation and his RelationXML
        inter->getNonSmoothLawPtr()->setNonSmoothLawXML(nslxml);

        // creation of the RelationXML
        static_cast<NewtonImpactFrictionNSLXML*>(nslxml)->updateNonSmoothLawXML(node, inter->getNonSmoothLawPtr());

        this->nSLawXML = nslxml;
      }
      else
      {
        XMLException::selfThrow("InteractionXML - loadInteraction error : undefined NonSmoothLaw type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
      }
    }
    else RuntimeException::selfThrow("InteractionXML::loadInteraction - There's no NonSmoothLaw in this Interaction, the XML platform can't be built");
  }
  else XMLException::selfThrow("InteractionXML - loadInteraction( Interaction* ) Error : no rootInteractionXMLNode is defined.");
  OUT("InteractionXML::loadInteraction( Interaction* )\n");
}

bool InteractionXML::hasAll() const
{
  if (SiconosDOMTreeTools::hasAttributeValue(dsConcernedNode, ALL_ATTRIBUTE))
    return SiconosDOMTreeTools::getBooleanAttributeValue(dsConcernedNode, ALL_ATTRIBUTE);
  else return false;
}

void InteractionXML::setAll(const bool& all)
{
  if (hasAll() == false)
  {
    if (all == true)
      xmlNewProp(dsConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str(), (xmlChar*)"true");
  }
  else if (all == false)
    xmlRemoveProp(xmlHasProp(dsConcernedNode, (xmlChar*)INTERACTION_DS_CONCERNED.c_str()));
}

