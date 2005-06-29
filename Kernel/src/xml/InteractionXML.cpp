#include "InteractionXML.h"

//---  Following includes to be suppressed thanks to factory (?) ---
// Relations
#include "LinearTIRXML.h"
#include "LagrangianLinearRXML.h"
#include "LagrangianNonLinearRXML.h"
// Nslaw
#include "ComplementarityConditionNSLXML.h"
#include "RelayNSLXML.h"
#include "NewtonImpactLawNSLXML.h"
#include "NewtonImpactFrictionNSLXML.h"
using namespace std;


InteractionXML::InteractionXML():
  rootInteractionXMLNode(NULL), idNode(NULL), nInterNode(NULL),
  statusNode(NULL),  yNode(NULL),   lambdaNode(NULL), isActiveNode(NULL),
  dsConcernedNode(NULL), dsListNode(NULL),  relationXML(NULL), nSLawXML(NULL),
  isRelationXMLAllocatedIn(false), isNsLawXMLAllocatedIn(false)
{}

InteractionXML::InteractionXML(xmlNode * interactionNode, vector<int> definedNumbers):
  rootInteractionXMLNode(interactionNode), idNode(NULL), nInterNode(NULL),
  statusNode(NULL),  yNode(NULL),   lambdaNode(NULL), isActiveNode(NULL),
  dsConcernedNode(NULL), dsListNode(NULL), relationXML(NULL), nSLawXML(NULL),
  isRelationXMLAllocatedIn(false), isNsLawXMLAllocatedIn(false)
{
  loadInteractionProperties(interactionNode, definedNumbers);
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

void InteractionXML::loadInteractionProperties(xmlNode * interactionNode, vector<int> definedNumbers)
{
  IN("InteractionXML::loadInteractionProperties\n ");
  xmlNode *node, *node2;
  string type;

  // id
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, ID_ATTRIBUTE)) != NULL)
    idNode = node;
  // nInter
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_NINTER)) != NULL)
    nInterNode = node;
  else
    XMLException::selfThrow("InteractionXML - loadInteractionProperties error : tag " + INTERACTION_NINTER + " not found.");
  // status
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_STATUS)) != NULL)
    statusNode = node;
  else
    XMLException::selfThrow("InteractionXML - loadInteractionProperties error : tag " + INTERACTION_STATUS + " not found.");
  // y
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_Y)) != NULL)
    yNode = node;
  // lambda
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_LAMBDA)) != NULL)
    lambdaNode = node;
  // dsConcerned
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_DS_CONCERNED)) != NULL)
  {
    dsConcernedNode = node;
    // Check if all ds are concerned or not
    if (! hasAll())
    {
      // Get the DSList node
      if (SiconosDOMTreeTools::findNodeChild(dsConcernedNode, INTERACTION_DS_LIST) != NULL)
        dsListNode = SiconosDOMTreeTools::findNodeChild(dsConcernedNode, INTERACTION_DS_LIST);
      else
        XMLException::selfThrow("tag DSlist not found.");
    }
  }
  else
    XMLException::selfThrow("InteractionXML - xml constructor, tag " + INTERACTION_DS_CONCERNED + " not found.");

  //  node = SiconosDOMTreeTools::findNodeChild((const xmlNode*)interactionNode, INTERACTION_CONTENT_TAG);
  if ((node = SiconosDOMTreeTools::findNodeChild((const xmlNode*)interactionNode, INTERACTION_CONTENT_TAG)) != NULL)
  {
    // the first child is the Relation
    if ((node2 = SiconosDOMTreeTools::findNodeChild(node)) != NULL)
    {
      // get Relation type
      type = (char*)node2->name;
      if (type == LINEAR_TIME_INVARIANT_RELATION_TAG)
        relationXML = new LinearTIRXML(node2);
      else if (type == LAGRANGIAN_NON_LINEAR_RELATION_TAG)
        relationXML = new LagrangianNonLinearRXML(node2);
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
  OUT("InteractionXML::loadInteractionProperties\n ");
}


void InteractionXML::loadInteractionConcernedDS(xmlNode * DSConcernedNode, vector<int> definedDSNumbers)
{
  xmlNode *DSnode;
  int number1, number2;
  //int size = SiconosDOMTreeTools::getIntegerAttributeValue(DSConcernedNode, INTERACTION_SIZE);
  unsigned int size = 0;
  unsigned int i = 0;
  vector<int> vtmp;

  if ((DSnode = SiconosDOMTreeTools::findNodeChild((const xmlNode*)DSConcernedNode, DYNAMICAL_SYSTEM_TAG)) == NULL)
    XMLException::selfThrow("InteractionXML - loadInteractionConcernedDS error : no " + DYNAMICAL_SYSTEM_TAG + " declared in " + INTERACTION_DS_CONCERNED + " tag.");

  size = SiconosDOMTreeTools::getNodeChildrenNumber(DSConcernedNode);
  while ((DSnode != NULL) && (i < size))
  {
    if ((number1 = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, NUMBER_ATTRIBUTE)) >= (number2 = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, INTERACTION_INTERACTWITHDS_NUMBER)))
      XMLException::selfThrow("InteractionXML - loadInteractionConcernedDS error : in tag " + INTERACTION_DS_CONCERNED + " a " + NUMBER_ATTRIBUTE + " is < to a " + INTERACTION_INTERACTWITHDS_NUMBER + " attribute.");

    //Verifying DS numbers couple exist
    unsigned int j = 0;
    while ((j < definedDSNumbers.size()) && (definedDSNumbers[j] != number1))
    {
      j++;
    }

    if (j == definedDSNumbers.size())
      XMLException::selfThrow("InteractionXML - loadInteractionConcernedDS error : in tag" + INTERACTION_DS_CONCERNED + " undefined DS number");
    j = 0;
    while ((j < definedDSNumbers.size()) && (definedDSNumbers[j] != number2))
    {
      j++;
    }

    if (i == definedDSNumbers.size())
      XMLException::selfThrow("InteractionXML - loadInteractionConcernedDS error : in tag" + INTERACTION_DS_CONCERNED + " undefined DS number");

    vtmp.clear();
    vtmp.push_back(number1);
    vtmp.push_back(number2);
    DSCouples.push_back(vtmp);

    if ((number1 = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, NUMBER_ATTRIBUTE)) >= (number2 = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, INTERACTION_INTERACTWITHDS_NUMBER)))
    {
      XMLException::selfThrow("InteractionXML - loadInteractionConcernedDS error : in a tag " + INTERACTION_DS_CONCERNED + " a " + NUMBER_ATTRIBUTE + " is < to a" + INTERACTION_INTERACTWITHDS_NUMBER + " attribute.");
    }

    DSnode = SiconosDOMTreeTools::findFollowNode(DSnode, DYNAMICAL_SYSTEM_TAG);
    i++;
  }

  if (i < size)
    XMLException::selfThrow("InteractionXML - loadInteractionConcernedDS error : wrong size attribute in the tag " + INTERACTION_DS_CONCERNED);
  OUT("InteractionXML::loadInteractionConcernedDS\n ");
}

void InteractionXML::updateInteractionXML(xmlNode* node, Interaction* inter)
{
  IN("InteractionXML::updateInteractionXML\n");
  this->rootInteractionXMLNode = node;
  this->loadInteraction(inter);
  OUT("InteractionXML::updateInteractionXML\n");
}

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
     * now, creation of the RelationXML object
     */
    xmlNode *InteractionContentNode;
    if (inter->getRelationPtr() != NULL)
    {
      type = inter->getRelationPtr()->getType();
      //node = xmlNewChild( this->rootInteractionXMLNode, NULL, (xmlChar*)INTERACTION_RELATION.c_str(), NULL );
      InteractionContentNode = xmlNewChild(this->rootInteractionXMLNode, NULL, (xmlChar*)INTERACTION_CONTENT_TAG.c_str(), NULL);
      if (type == LAGRANGIANLINEARRELATION)
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
      else if (type == LAGRANGIANNONLINEARRELATION)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_LNL.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)LAGRANGIAN_NON_LINEAR_RELATION_TAG.c_str(), NULL);
        newRelationXml = new LagrangianNonLinearRXML();

        // linkage between the Relation and his RelationXML
        inter->getRelationPtr()->setRelationXML(newRelationXml);

        // creation of the RelationXML
        static_cast<LagrangianNonLinearRXML*>(newRelationXml)->updateRelationXML(node, inter->getRelationPtr());

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
      else if (type == NEWTONIMPACTLAWNSLAW)
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

void InteractionXML::setDSConcerned(vector<DynamicalSystem*> dsConcerned)
{
  if (this->dsConcernedNode == NULL)
  {
    unsigned int i;
    unsigned int number1, number2;
    vector<int> vtmp;
    xmlNode* node;

    // conversion of the vector<DynamicalSystem*> in vector< vector<int> >
    this->DSCouples.clear();
    for (i = 0; i < dsConcerned.size(); i = i + 2)
    {
      number1 = dsConcerned[i]->getNumber();
      if ((i + 1) < (dsConcerned.size()))
        number2 = dsConcerned[i + 1]->getNumber();
      else cout << "Warning: uneven number of dynamical systems in the interaction" << endl;
      if (number1 > number2)
      {
        vtmp.clear();
        vtmp.push_back(number2);
        vtmp.push_back(number1);
      }
      else if (number1 < number2)
      {
        vtmp.clear();
        vtmp.push_back(number1);
        vtmp.push_back(number2);
      }
      else
        XMLException::selfThrow("InteractionXML - loadInteraction ERROR : This interaction is about the same dynamical system.");

      this->DSCouples.push_back(vtmp);
    }

    //save in the DOM tree
    char num[32];

    /*
     * creation of the DS_Concerned node
     */
    this->dsConcernedNode = xmlNewChild(this->rootInteractionXMLNode, NULL, (xmlChar*)INTERACTION_DS_CONCERNED.c_str(), NULL);
    //sprintf(num, "%i", this->DSCouples.size());
    //xmlNewProp(this->dsConcernedNode, (xmlChar*)INTERACTION_SIZE.c_str(), (xmlChar*)num);

    for (i = 0; i < this->DSCouples.size(); i++)
    {
      node = xmlNewChild(this->dsConcernedNode, NULL, (xmlChar*)DYNAMICAL_SYSTEM_TAG.c_str(), NULL);
      sprintf(num, "%i", this->DSCouples[i][0]);
      xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
      sprintf(num, "%i", this->DSCouples[i][1]);
      xmlNewProp(node, (xmlChar*)INTERACTION_INTERACTWITHDS_NUMBER.c_str(), (xmlChar*)num);
    }
  }
  else
  {
    /* \todo : when interaction have been given in the XML input file and that the user has added new ones */
  }
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

