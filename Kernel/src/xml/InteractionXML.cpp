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
#include "InteractionXML.hpp"

//---  Following includes to be suppressed thanks to factory (?) ---
// Relations
#include "LinearRXML.hpp"
// Nslaw
#include "ComplementarityConditionNSLXML.hpp"
#include "RelayNSLXML.hpp"
#include "NewtonImpactNSLXML.hpp"
#include "NewtonImpactFrictionNSLXML.hpp"

#include "NonSmoothLaw.hpp"
using namespace std;


InteractionXML::InteractionXML():
  rootNode(NULL), sizeNode(NULL), yNode(NULL), lambdaNode(NULL), DSConcernedNode(NULL),  relationNode(NULL), nsLawNode(NULL)
{}

InteractionXML::InteractionXML(xmlNodePtr  interactionNode):
  rootNode(interactionNode), sizeNode(NULL), yNode(NULL), lambdaNode(NULL), DSConcernedNode(NULL),  relationNode(NULL), nsLawNode(NULL)
{
  xmlNodePtr node, node2;
  string type;
  string subType;
  // size (required)
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, "size")))
    sizeNode = node;
  else
    XMLException::selfThrow("InteractionXML - xml constructor: tag size not found.");

  // y (optional)
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, "y")))
    yNode = node;
  // lambda (optional)
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, "lambda")))
    lambdaNode = node;
  // get DSConcerned node (required)
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, DS_CONCERNED)))
    DSConcernedNode = node;
  else
    XMLException::selfThrow("InteractionXML - xml constructor, tag " + DS_CONCERNED + " not found.");

  // Relation and Non Smooth Law (required)
  if ((node = SiconosDOMTreeTools::findNodeChild((const xmlNodePtr)interactionNode, INTERACTION_CONTENT_TAG)))
  {
    // the first child is the Relation
    if ((node2 = SiconosDOMTreeTools::findNodeChild(node)))
    {
      relationNode = node2;
      // get Relation sub-type
      subType = SiconosDOMTreeTools::getStringAttributeValue(node2, "type");

      // Linear case
      if (subType == "Linear" || subType == "LinearTI")
        relationXML.reset(new LinearRXML(node2));
      // General non linear case
      else
        relationXML.reset(new RelationXML(node2));
    }
    else
      XMLException::selfThrow("InteractionXML - constructor: relation not found.");

    // the second child is the NonSmoothLaw
    nsLawNode = SiconosDOMTreeTools::findFollowNode(node2);
    if (nsLawNode)
    {
      // Non smooth law type
      type = (char*)nsLawNode->name;
      if (type == COMPLEMENTARITY_CONDITION_NSLAW_TAG)
        nSLawXML.reset(new ComplementarityConditionNSLXML(nsLawNode));
      else if (type == RELAY_NSLAW_TAG)
        nSLawXML.reset(new RelayNSLXML(nsLawNode));
      else if (type == NEWTON_IMPACT_NSLAW_TAG)
        nSLawXML.reset(new NewtonImpactNSLXML(nsLawNode));
      else if (type == NEWTON_IMPACT_FRICTION_NSLAW_TAG)
        nSLawXML.reset(new NewtonImpactFrictionNSLXML(nsLawNode));
      else
        XMLException::selfThrow("InteractionXML, constructor(xml) undefined non smooth law type : " + type);
    }
    else
      XMLException::selfThrow("InteractionXML - constructor error: NonSmoothLaw not found.");
  }
}

InteractionXML::~InteractionXML()
{
}

void InteractionXML::setSize(const unsigned int newSize)
{
  if (!hasSize())
    sizeNode = SiconosDOMTreeTools::createIntegerNode(rootNode, "size", newSize);
  else SiconosDOMTreeTools::setIntegerContentValue(sizeNode, newSize);
}

void InteractionXML::setY(const SiconosVector& v)
{
  if (!hasY())
    yNode = SiconosDOMTreeTools::createVectorNode(rootNode, "y", v);
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(yNode, v);
}

void InteractionXML::setLambda(const SiconosVector& v)
{
  if (!hasLambda())
    lambdaNode = SiconosDOMTreeTools::createVectorNode(rootNode, "lambda", v);
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(lambdaNode, v);
}

void InteractionXML::getDSNumbers(vector<int>& dsNumbers)
{
  if (!hasAllDS())
    SiconosDOMTreeTools::getVector(DSConcernedNode, dsNumbers);
  else
    XMLException::selfThrow("InteractionXML::getDSNumbers - No list of ds, all parameter = true.");
}


void InteractionXML::updateInteractionXML(xmlNodePtr  node, SP::Interaction inter)
{
  rootNode = node;
  loadInteraction(inter);
}

// warning: this function has not been reviewed for multiple DS loading !!
void InteractionXML::loadInteraction(SP::Interaction inter)
{
  XMLException::selfThrow("InteractionXML::loadInteraction - not implemented.");

  //   string type;
  //   string tmp;
  //   xmlNodePtr  node;
  //   RelationXML* newRelationXml;
  //   NonSmoothLawXML* nslxml;

  //   if( this->rootNode != NULL )
  //     {
  //       /*
  //        * Creation of the RelationXML object
  //        */
  //       xmlNodePtr InteractionContentNode;
  //       if( inter->relation() != NULL )
  //  {
  //    type = inter->relation()->getType();
  //    InteractionContentNode = xmlNewChild( this->rootNode, NULL, (xmlChar*)INTERACTION_CONTENT_TAG.c_str(), NULL );

  //    if (type == RELATION)
  //      {
  //        node = xmlNewChild( InteractionContentNode, NULL, (xmlChar*)RELATION_TAG.c_str(), NULL );
  //        relationXML = new RelationXML(node);
  //        // linkage between the Relation and his RelationXML
  //        inter->relation()->setRelationXML( relationXML );
  //      }
  //    else if (type == LAGRANGIANLINEARRELATION)
  //      {
  //        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_LL.c_str() );
  //        node = xmlNewChild( InteractionContentNode, NULL, (xmlChar*)LAGRANGIAN_LINEAR_RELATION_TAG.c_str(), NULL );
  //        newRelationXml = new LagrangianLinearRXML();

  //        // linkage between the Relation and his RelationXML
  //        inter->relation()->setRelationXML( newRelationXml );

  //        // creation of the RelationXML
  //        static_cast<LagrangianLinearRXML*>(newRelationXml)->updateRelationXML( node, inter->relation() );

  //        relationXML = newRelationXml;
  //      }
  //    else if (type == LAGRANGIANRELATION)
  //      {
  //        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_LNL.c_str() );
  //        node = xmlNewChild( InteractionContentNode, NULL, (xmlChar*)LAGRANGIAN_RELATION_TAG.c_str(), NULL );
  //        newRelationXml = new LagrangianRXML();

  //        // linkage between the Relation and his RelationXML
  //        inter->relation()->setRelationXML( newRelationXml );

  //        // creation of the RelationXML
  //        static_cast<LagrangianRXML*>(newRelationXml)->updateRelationXML( node, inter->relation() );

  //        relationXML = newRelationXml;
  //      }
  //    else if (type == LINEARTIRELATION)
  //      {
  //        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_LTI.c_str() );
  //        node = xmlNewChild( InteractionContentNode, NULL, (xmlChar*)LINEAR_TIME_INVARIANT_RELATION_TAG.c_str(), NULL );
  //        newRelationXml = new LinearTIRXML();

  //        // linkage between the Relation and his RelationXML
  //        inter->relation()->setRelationXML( newRelationXml );

  //        // creation of the RelationXML
  //        static_cast<LinearTIRXML*>(newRelationXml)->updateRelationXML( node, inter->relation() );

  //        relationXML = newRelationXml;
  //      }
  //    else
  //      {
  //        XMLException::selfThrow("InteractionXML - loadInteraction error : undefined Relation type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
  //      }
  //  }
  //       else RuntimeException::selfThrow("InteractionXML::loadInteraction - There's no Relation in this Interaction, the XML platform can't be built");


  //       /*
  //        * now, creation of the NonSmoothLawXML object
  //        */
  //       if( inter->nonSmoothLaw() != NULL )
  //  {
  //    type = inter->nonSmoothLaw()->getType();
  //    //node = xmlNewChild( this->rootNode, NULL, (xmlChar*)INTERACTION_NS_LAW.c_str(), NULL );
  //    if (type == COMPLEMENTARITYCONDITIONNSLAW)
  //      {
  //        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_COMPLEMENTARITYCONDITIONNSLAW.c_str() );
  //        node = xmlNewChild( InteractionContentNode, NULL, (xmlChar*)COMPLEMENTARITY_CONDITION_NSLAW_TAG.c_str(), NULL );
  //        nslxml = new ComplementarityConditionNSLXML();

  //        // linkage between the Relation and his RelationXML
  //        inter->nonSmoothLaw()->setNonSmoothLawXML( nslxml );

  //        // creation of the RelationXML
  //        static_cast<ComplementarityConditionNSLXML*>(nslxml)->updateNonSmoothLawXML( node, inter->nonSmoothLaw() );

  //        this->nSLawXML = nslxml;
  //      }
  //    else if (type == RELAYNSLAW)
  //      {
  //        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_RELAYNSLAW.c_str() );
  //        node = xmlNewChild( InteractionContentNode, NULL, (xmlChar*)RELAY_NSLAW_TAG.c_str(), NULL );
  //        nslxml = new RelayNSLXML();

  //        // linkage between the Relation and his RelationXML
  //        inter->nonSmoothLaw()->setNonSmoothLawXML( nslxml );

  //        // creation of the RelationXML
  //        static_cast<RelayNSLXML*>(nslxml)->updateNonSmoothLawXML( node, inter->nonSmoothLaw() );

  //        this->nSLawXML = nslxml;
  //      }
  //    else if (type == NEWTONIMPACTNSLAW)
  //      {
  //        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_NEWTONIMPACTLAWNSLAW.c_str() );
  //        node = xmlNewChild( InteractionContentNode, NULL, (xmlChar*)NEWTON_IMPACT_LAW_NSLAW_TAG.c_str(), NULL );
  //        nslxml = new NewtonImpactNSLXML();

  //        // linkage between the Relation and his RelationXML
  //        inter->nonSmoothLaw()->setNonSmoothLawXML( nslxml );

  //        // creation of the RelationXML
  //        static_cast<NewtonImpactNSLXML*>(nslxml)->updateNonSmoothLawXML( node, inter->nonSmoothLaw() );

  //        this->nSLawXML = nslxml;
  //      }
  //    else if (type == NEWTONIMPACTFRICTIONNSLAW)
  //      {
  //        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_NEWTONIMPACTLAWNSLAW.c_str() );
  //        node = xmlNewChild( InteractionContentNode, NULL, (xmlChar*)NEWTON_IMPACT_FRICTION_NSLAW_TAG.c_str(), NULL );
  //        nslxml = new NewtonImpactFrictionNSLXML();

  //        // linkage between the Relation and his RelationXML
  //        inter->nonSmoothLaw()->setNonSmoothLawXML( nslxml );

  //        // creation of the RelationXML
  //        static_cast<NewtonImpactFrictionNSLXML*>(nslxml)->updateNonSmoothLawXML( node, inter->nonSmoothLaw() );

  //        this->nSLawXML = nslxml;
  //      }
  //    else
  //      {
  //        XMLException::selfThrow("InteractionXML - loadInteraction error : undefined NonSmoothLaw type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
  //      }
  //  }
  //       else RuntimeException::selfThrow("InteractionXML::loadInteraction - There's no NonSmoothLaw in this Interaction, the XML platform can't be built");
  //     }
  //   else XMLException::selfThrow("InteractionXML - loadInteraction( SP::Interaction ) Error : no rootNode is defined.");
}

