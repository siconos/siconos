//$Id: InteractionXML.cpp,v 1.43 2005/03/22 15:55:05 jbarbier Exp $
#include "InteractionXML.h"

#include "LinearTIRXML.h"
#include "LagrangianLinearRXML.h"
#include "LagrangianNonLinearRXML.h"
#include "ComplementarityConditionNSLXML.h"
#include "RelayNSLXML.h"
#include "NewtonImpactLawNSLXML.h"
#include "NewtonImpactFrictionNSLXML.h"

#include "check.h"


InteractionXML::InteractionXML()
{
  this->idNode = NULL;
  this->nInterNode = NULL;
  this->statusNode = NULL;
  this->yNode = NULL;
  this->lambdaNode = NULL;
  this->isActiveNode = NULL;
  this->dsConcernedNode = NULL;
}

InteractionXML::InteractionXML(xmlNode * interactionNode, vector<int> definedNumbers)
{
  this->idNode = NULL;
  this->nInterNode = NULL;
  this->statusNode = NULL;
  this->yNode = NULL;
  this->lambdaNode = NULL;
  this->isActiveNode = NULL;
  this->dsConcernedNode = NULL;
  this->rootInteractionXMLNode = interactionNode;

  loadInteractionProperties(interactionNode, definedNumbers);
}

InteractionXML::~InteractionXML()
{
  //  if(this->relationXML != NULL) delete this->relationXML;
  //  if(this->nSLawXML != NULL) delete this->nSLawXML;
}

void InteractionXML::loadInteractionProperties(xmlNode * interactionNode, vector<int> definedNumbers)
{
  IN("InteractionXML::loadInteractionProperties\n ");
  xmlNode *node;
  string type;

  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, ID_ATTRIBUTE)) != NULL)
  {
    this->idNode = node;
  }
  else
  {
    //XMLException::selfThrow("InteractionXML - loadInteractionProperties error : tag " + INTERACTION_ID + " not found.");
    cout << "InteractionXML - loadInteractionProperties WARNING : tag " << ID_ATTRIBUTE << " not found, this attribute is optional." << endl;
    this->idNode = NULL;
  }
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_NINTER)) != NULL)
  {
    this->nInterNode = node;
  }
  else
  {
    XMLException::selfThrow("InteractionXML - loadInteractionProperties error : tag " + INTERACTION_NINTER + " not found.");
  }
  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_STATUS)) != NULL)
  {
    this->statusNode = node;
  }
  else
  {
    XMLException::selfThrow("InteractionXML - loadInteractionProperties error : tag " + INTERACTION_STATUS + " not found.");
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_Y)) != NULL)
  {
    this->yNode = node;
  }
  else
  {
    cout << "InteractionXML - loadInteractionProperties WARNING : tag " << INTERACTION_Y << " not found, this attribute is optional." << endl;
    this->yNode = NULL;
  }


  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_LAMBDA)) != NULL)
  {
    this->lambdaNode = node;
  }
  else
  {
    cout << "InteractionXML - loadInteractionProperties WARNING : tag " << INTERACTION_LAMBDA << " not found, this attribute is optional." << endl;
    this->lambdaNode = NULL;
  }
  //    if ((node=SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_ISACTIVE)) !=NULL)
  //    {
  //    this->isActiveNode=node;
  //    }
  //    else
  //    {
  //    //XMLException::selfThrow("InteractionXML - loadInteractionProperties error : tag " + INTERACTION_ISACTIVE + " not found.");
  //    this->isActiveNode=NULL;
  //    }

  if ((node = SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_DS_CONCERNED)) != NULL)
  {
    this->dsConcernedNode = node;
    loadInteractionConcernedDS(node, definedNumbers);
  }
  else
  {
    XMLException::selfThrow("InteractionXML - loadInteractionProperties error : tag " + INTERACTION_DS_CONCERNED + " not found.");
  }

  node = SiconosDOMTreeTools::findNodeChild((const xmlNode*)interactionNode, INTERACTION_CONTENT_TAG);
  //if ((node=SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_RELATION)) !=NULL)
  if (node != NULL)
  {
    /*
     * the first child is the Relation
     */
    node = SiconosDOMTreeTools::findNodeChild(node);
    if (node != NULL)
    {
      //string type = SiconosDOMTreeTools::getStringAttributeValue(node,INTERACTION_TYPE); //type of Relation
      type = (char*)node->name;
      if (type == LINEAR_TIME_INVARIANT_RELATION_TAG)
      {
        this->relationXML = new LinearTIRXML(node);
      }
      else if (type == LAGRANGIAN_NON_LINEAR_RELATION_TAG)
      {
        this->relationXML = new LagrangianNonLinearRXML(node);
      }
      else if (type == LAGRANGIAN_LINEAR_RELATION_TAG)
      {
        this->relationXML = new LagrangianLinearRXML(node);
      }
      else
      {
        XMLException::selfThrow("InteractionXML : undefined Relation type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?)");
      }
    }
    else
    {
      XMLException::selfThrow("InteractionXML - loadInteractionProperties error : Relation not found.");
    }

    /*
       * the second child is the NonSmoothLaw
       */
    node = SiconosDOMTreeTools::findFollowNode(node);
    //if ((node=SiconosDOMTreeTools::findNodeChild(interactionNode, INTERACTION_NS_LAW)) !=NULL)
    if (node != NULL)
    {
      //string type = SiconosDOMTreeTools::getStringAttributeValue(node,INTERACTION_TYPE); //type of Relation
      type = (char*)node->name;
      if (type == COMPLEMENTARITY_CONDITION_NSLAW_TAG)
      {
        this->nSLawXML = new ComplementarityConditionNSLXML(node);
      }
      else if (type == RELAY_NSLAW_TAG)
      {
        this->nSLawXML = new RelayNSLXML(node);
      }
      else if (type == NEWTON_IMPACT_LAW_NSLAW_TAG)
      {
        this->nSLawXML = new NewtonImpactLawNSLXML(node);
      }
      else if (type == NEWTON_IMPACT_FRICTION_NSLAW_TAG)
      {
        this->nSLawXML = new NewtonImpactFrictionNSLXML(node);
      }
      else
      {
        XMLException::selfThrow("InteractionXML : undefined NonSmoothLaw type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?)");
      }
    }
    else
    {
      XMLException::selfThrow("InteractionXML - loadInteractionProperties error : NonSmoothLaw not found.");
    }
  }

  OUT("InteractionXML::loadInteractionProperties\n ");

}


void InteractionXML::loadInteractionConcernedDS(xmlNode * DSConcernedNode, vector<int> definedDSNumbers)
{
  IN("InteractionXML::loadInteractionConcernedDS\n ");

  xmlNode *DSnode;
  int number1, number2;
  //int size = SiconosDOMTreeTools::getIntegerAttributeValue(DSConcernedNode, INTERACTION_SIZE);
  int size = 0;
  int i = 0;
  vector<int> vtmp;

  if ((DSnode = SiconosDOMTreeTools::findNodeChild((const xmlNode*)DSConcernedNode, DYNAMICAL_SYSTEM_TAG)) == NULL)
  {
    XMLException::selfThrow("InteractionXML - loadInteractionConcernedDS error : at least one couple of " + DYNAMICAL_SYSTEM_TAG + " must be declared in " + INTERACTION_DS_CONCERNED + " tag.");
  }

  size = SiconosDOMTreeTools::getNodeChildrenNumber(DSConcernedNode);
  //  cout<<"number of children = "<<size<<endl;
  //  cout<<"# Press <<ENTER>> !!"<<endl;
  //  getchar();

  while ((DSnode != NULL) && (i < size))
  {
    if ((number1 = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, NUMBER_ATTRIBUTE)) >= (number2 = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, INTERACTION_INTERACTWITHDS_NUMBER)))
    {
      XMLException::selfThrow("InteractionXML - loadInteractionConcernedDS error : in a tag " + INTERACTION_DS_CONCERNED + " a " + NUMBER_ATTRIBUTE + " is < to a " + INTERACTION_INTERACTWITHDS_NUMBER + " attribute.");
    }

    //  cout<<"number1 = "<<number1<<"  --   number2 = "<<number2<<endl;
    //  cout<<"# Press <<ENTER>> !!"<<endl;
    //  getchar();


    //Verifying DS numbers couple exist
    int j = 0;
    while ((j < definedDSNumbers.size()) && (definedDSNumbers[j] != number1))
    {
      j++;
    }

    if (j == definedDSNumbers.size())
    {
      char errorMsg[1024];
      sprintf(errorMsg, "InteractionXML - loadInteractionConcernedDS error : in a tag %s you define couple of DS with a DS number who doesn't exist : %d.", INTERACTION_DS_CONCERNED.c_str(), number1);
      XMLException::selfThrow(errorMsg);
    }
    j = 0;
    while ((j < definedDSNumbers.size()) && (definedDSNumbers[j] != number2))
    {
      j++;
    }

    if (i == definedDSNumbers.size())
    {
      char errorMsg[1024];
      sprintf(errorMsg, "InteractionXML - loadInteractionConcernedDS error : in a tag %s you define couple of DS with a DS number who doesn't exist : %d.", INTERACTION_DS_CONCERNED.c_str(), number2);
      XMLException::selfThrow(errorMsg);
    }


    vtmp.clear();
    vtmp.push_back(number1);
    vtmp.push_back(number2);
    DSCouples.push_back(vtmp);

    if ((number1 = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, NUMBER_ATTRIBUTE)) >= (number2 = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, INTERACTION_INTERACTWITHDS_NUMBER)))
    {
      XMLException::selfThrow("InteractionXML - loadInteractionConcernedDS error : in a tag " + INTERACTION_DS_CONCERNED + " a " + NUMBER_ATTRIBUTE + " is < to a" + INTERACTION_INTERACTWITHDS_NUMBER + " attribute.");
    }

    //    //Verifying DS numbers couple exist
    //    j=0;
    //    while ((j<definedDSNumbers.size()) && (definedDSNumbers[j]!=number1))
    //    {
    //      j++;
    //    }
    //    if (j==definedDSNumbers.size())
    //    {
    //      char errorMsg[1024];
    //      sprintf(errorMsg, "InteractionXML - loadInteractionConcernedDS error : in a tag %s you define couple of DS with a DS number who doesn't exist : %d.", INTERACTION_DS_CONCERNED.c_str(), number1);
    //      XMLException::selfThrow(errorMsg);
    //    }
    //
    //    j=0;
    //    while ((j<definedDSNumbers.size()) && (definedDSNumbers[j]!=number2))
    //    {
    //      j++;
    //    }
    //    if (i==definedDSNumbers.size())
    //    {
    //      char errorMsg[1024];
    //      sprintf(errorMsg, "InteractionXML - loadInteractionConcernedDS error : in a tag %s you define couple of DS with a DS number who doesn't exist : %d.", INTERACTION_DS_CONCERNED.c_str(), number2);
    //      XMLException::selfThrow(errorMsg);
    //    }
    //
    //      this->DSCouples[i][0]=number1;
    //    this->DSCouples[i][1]=number2;

    DSnode = SiconosDOMTreeTools::findFollowNode(DSnode, DYNAMICAL_SYSTEM_TAG);

    i++;
  }

  if (i < size)
  {
    XMLException::selfThrow("InteractionXML - loadInteractionConcernedDS error : the size attribute given in the tag " + INTERACTION_DS_CONCERNED + " is not correct.");
  }
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
  RelationXML* relationxml;
  NonSmoothLawXML* nslxml;
  int number, i;
  char num[32];

  if (this->rootInteractionXMLNode != NULL)
  {
    /*
     * now, creation of the RelationXML object
     */
    xmlNode *InteractionContentNode;
    if (inter->getRelation() != NULL)
    {
      type = inter->getRelation()->getType();
      //node = xmlNewChild( this->rootInteractionXMLNode, NULL, (xmlChar*)INTERACTION_RELATION.c_str(), NULL );
      InteractionContentNode = xmlNewChild(this->rootInteractionXMLNode, NULL, (xmlChar*)INTERACTION_CONTENT_TAG.c_str(), NULL);
      if (type == LAGRANGIANLINEARRELATION)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_LL.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)LAGRANGIAN_LINEAR_RELATION_TAG.c_str(), NULL);
        relationxml = new LagrangianLinearRXML();

        // linkage between the Relation and his RelationXML
        inter->getRelation()->setRelationXML(relationxml);

        // creation of the RelationXML
        static_cast<LagrangianLinearRXML*>(relationxml)->updateRelationXML(node, inter->getRelation());

        this->relationXML = relationxml;
      }
      else if (type == LAGRANGIANNONLINEARRELATION)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_LNL.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)LAGRANGIAN_NON_LINEAR_RELATION_TAG.c_str(), NULL);
        relationxml = new LagrangianNonLinearRXML();

        // linkage between the Relation and his RelationXML
        inter->getRelation()->setRelationXML(relationxml);

        // creation of the RelationXML
        static_cast<LagrangianNonLinearRXML*>(relationxml)->updateRelationXML(node, inter->getRelation());

        this->relationXML = relationxml;
      }
      else if (type == LINEARTIRELATION)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_LTI.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)LINEAR_TIME_INVARIANT_RELATION_TAG.c_str(), NULL);
        relationxml = new LinearTIRXML();

        // linkage between the Relation and his RelationXML
        inter->getRelation()->setRelationXML(relationxml);

        // creation of the RelationXML
        static_cast<LinearTIRXML*>(relationxml)->updateRelationXML(node, inter->getRelation());

        this->relationXML = relationxml;
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
    if (inter->getNonSmoothLaw() != NULL)
    {
      type = inter->getNonSmoothLaw()->getType();
      //node = xmlNewChild( this->rootInteractionXMLNode, NULL, (xmlChar*)INTERACTION_NS_LAW.c_str(), NULL );
      if (type == COMPLEMENTARITYCONDITIONNSLAW)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_COMPLEMENTARITYCONDITIONNSLAW.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)COMPLEMENTARITY_CONDITION_NSLAW_TAG.c_str(), NULL);
        nslxml = new ComplementarityConditionNSLXML();

        // linkage between the Relation and his RelationXML
        inter->getNonSmoothLaw()->setNonSmoothLawXML(nslxml);

        // creation of the RelationXML
        static_cast<ComplementarityConditionNSLXML*>(nslxml)->updateNonSmoothLawXML(node, inter->getNonSmoothLaw());

        this->nSLawXML = nslxml;
      }
      else if (type == RELAYNSLAW)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_RELAYNSLAW.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)RELAY_NSLAW_TAG.c_str(), NULL);
        nslxml = new RelayNSLXML();

        // linkage between the Relation and his RelationXML
        inter->getNonSmoothLaw()->setNonSmoothLawXML(nslxml);

        // creation of the RelationXML
        static_cast<RelayNSLXML*>(nslxml)->updateNonSmoothLawXML(node, inter->getNonSmoothLaw());

        this->nSLawXML = nslxml;
      }
      else if (type == NEWTONIMPACTLAWNSLAW)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_NEWTONIMPACTLAWNSLAW.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)NEWTON_IMPACT_LAW_NSLAW_TAG.c_str(), NULL);
        nslxml = new NewtonImpactLawNSLXML();

        // linkage between the Relation and his RelationXML
        inter->getNonSmoothLaw()->setNonSmoothLawXML(nslxml);

        // creation of the RelationXML
        static_cast<NewtonImpactLawNSLXML*>(nslxml)->updateNonSmoothLawXML(node, inter->getNonSmoothLaw());

        this->nSLawXML = nslxml;
      }
      else if (type == NEWTONIMPACTFRICTIONNSLAW)
      {
        //xmlNewProp( node, (xmlChar*)INTERACTION_TYPE.c_str(), (xmlChar*)INTERACTION_NEWTONIMPACTLAWNSLAW.c_str() );
        node = xmlNewChild(InteractionContentNode, NULL, (xmlChar*)NEWTON_IMPACT_FRICTION_NSLAW_TAG.c_str(), NULL);
        nslxml = new NewtonImpactFrictionNSLXML();

        // linkage between the Relation and his RelationXML
        inter->getNonSmoothLaw()->setNonSmoothLawXML(nslxml);

        // creation of the RelationXML
        static_cast<NewtonImpactFrictionNSLXML*>(nslxml)->updateNonSmoothLawXML(node, inter->getNonSmoothLaw());

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
    int i;
    int number1, number2;
    vector<int> vtmp;
    xmlNode* node;

    // conversion of the vector<DynamicalSystem*> in vector< vector<int> >
    this->DSCouples.clear();
    for (i = 0; i < dsConcerned.size(); i = i + 2)
    {
      number1 = dsConcerned[i]->getNumber();
      if ((i + 1) < (dsConcerned.size()))
        number2 = dsConcerned[i + 1]->getNumber();
      else cout << "!!! /!\ the number of dynamical system of the interaction is unpair, a couple of dynamical system in interaction is not complete." << endl;

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

//$Log: InteractionXML.cpp,v $
//Revision 1.43  2005/03/22 15:55:05  jbarbier
//- class NewtonImpactFriction non smooth law added to the kernel
//
//- xml schema modified for this new class
//- xml schema modified to accept a "joker" for further use of a LMGC90 mechanical plugin
//
//- new test added for the loading/saving of a NewtonImpactFrictionNSL
//
//Revision 1.42  2005/03/10 12:55:21  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.41  2005/03/09 15:30:36  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.40  2005/03/08 12:41:38  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.39  2005/03/07 13:17:21  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.38  2005/03/04 15:35:27  jbarbier
//- README files added for some samples
//
//- beginning of the refactoring of XML module constants
//
//Revision 1.37  2005/01/13 14:14:39  jbarbier
//- correction in the XML output about size attribute in tags DS_Concerned and Interactoin _Concerned
//
//- modifications in creation of XML objects when saving data with partial XML input file
//
//Revision 1.36  2005/01/10 17:06:37  jbarbier
//- attribute "size" is now unused in the code
//
//- xml schema v1.2 is in progress
//
//Revision 1.35  2004/12/08 12:49:38  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.34  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.33  2004/09/14 13:49:55  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.32  2004/09/10 08:04:49  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.31  2004/07/29 14:04:00  jbarbier
//- new test on SiconosMemoryXML
//
//- last functions hasAttribute() in the XML part added
//
//Revision 1.30  2004/07/06 14:54:49  acary
//Renaming NSLaw into NonSmoothLaw
//Renaming RelayNSLaw into RelayNSL
//Renaming CCNSLaw into ComplementarityConditionNSL
//Renaming NewtonImpactLaw into NewtonImpactLawNSL
//
//Revision 1.29  2004/07/06 11:48:09  acary
//Renaming LTIRelation into LinearTIR
//
//Revision 1.28  2004/07/06 10:25:20  acary
//Renaming LNLRelation into LagrangianNonLinearR
//
//Revision 1.27  2004/07/06 10:00:01  acary
//Change naming of LagrangianLineaR
//
//Revision 1.26  2004/07/06 08:09:10  acary
//Renaming Class LLRelation into LagrangianLinearR
//
//Revision 1.25  2004/06/30 09:44:35  acary
//Added NewtonImpactLawNSL
//