//$Id: NSDSXML.cpp,v 1.47 2005/03/23 15:03:55 jbarbier Exp $

#include "NSDSXML.h"

#include "LagrangianNLDSXML.h"
#include "LagrangianTIDSXML.h"
#include "LinearSystemDSXML.h"

#include "LinearECXML.h"
#include "LinearTIECXML.h"
#include "LagrangianECXML.h"

#include "check.h"


NSDSXML::NSDSXML()
{
  this->NSDSNode = NULL;
  this->DSXMLMap.clear();
}

NSDSXML::NSDSXML(xmlNode * rootNSDSNode)
{
  IN("NSDSXML::NSDSXML(xmlNode * rootNSDSNode)\n");
  if (rootNSDSNode != NULL)
  {
    this->NSDSNode = rootNSDSNode;
    this->DSXMLMap.clear();
    loadNonSmoothDynamicalSystem();
  }
  OUT("NSDSXML::NSDSXML(xmlNode * rootNSDSNode)\n");
}

NSDSXML::~NSDSXML()
{
  if (this->DSXMLMap.size() > 0)
  {
    for (int i = 0; i < this->DSXMLMap.size(); i++)
    {
      delete this->DSXMLMap[i];
    }
    this->DSXMLMap.clear();
  }

  if (this->interactionXMLMap.size() > 0)
  {
    for (int i = 0; i < this->interactionXMLMap.size(); i++)
    {
      delete this->interactionXMLMap[i];
    }
    this->interactionXMLMap.clear();
  }
}


/* get the DS of number "number"*/
DSXML* NSDSXML::getDSXML(int number)
{
  map<int, DSXML*>::iterator it;

  it = this->DSXMLMap.find(number);
  if (it == this->DSXMLMap.end())
  {
    cout << "NSDSXML::getDSXML - Error : the DSXML number " << number << " does not exist!" << endl;
    return NULL;
  }
  return this->DSXMLMap[number];
}


InteractionXML* NSDSXML::getInteractionXML(int number)
{
  map<int, InteractionXML*>::iterator it;

  it = this->interactionXMLMap.find(number);
  if (it == this->interactionXMLMap.end())
  {
    cout << "NSDSXML::getInteractionXML - Error : the InteractionXML number " << number << " does not exist!" << endl;
    return NULL;
  }
  return this->interactionXMLMap[number];
}

EqualityConstraintXML* NSDSXML::getEqualityConstraintXML(int number)
{
  map<int, EqualityConstraintXML*>::iterator it;

  it = this->equalityConstraintXMLMap.find(number);
  if (it == this->equalityConstraintXMLMap.end())
  {
    cout << "NSDSXML::getEqualityConstraintXML - Error : the EqualityConstraintXML number " << number << " does not exist!" << endl;
    return NULL;
  }
  return this->equalityConstraintXMLMap[number];
}


void NSDSXML::loadNonSmoothDynamicalSystem()
{
  xmlNode *node;
  if ((node = SiconosDOMTreeTools::findNodeChild(NSDSNode, LMGC90_NSDS_TAG)) == NULL)
  {
    if ((node = SiconosDOMTreeTools::findNodeChild(NSDSNode, DYNAMICAL_SYSTEM_DEFINITION_TAG)) != NULL)
      this->loadDSXML(node);
    else
      XMLException::selfThrow("NSDSXML - loadNSDS error : tag " + DYNAMICAL_SYSTEM_DEFINITION_TAG + " not found.");

    if ((node = SiconosDOMTreeTools::findNodeChild(NSDSNode, INTERACTION_DEFINITION_TAG)) != NULL)
    {
      this->loadInteractionXML(node);
    }
    else
      //XMLException::selfThrow("NSDSXML - loadNSDS error : tag " + NSDS_INTERACTION_DEFINITION + " not found.");
      cout << "NSDSXML - loadNSDS WARNING : tag " << INTERACTION_DEFINITION_TAG << " not found,\nDefining interactions is optional." << endl;

    if ((node = SiconosDOMTreeTools::findNodeChild(NSDSNode, EQUALITYCONSTRAINT_DEFINITION_TAG)) != NULL)
    {
      this->loadEqualityConstraintXML(node);
    }
    else
      cout << "NSDSXML - loadNSDS WARNING : tag " << EQUALITYCONSTRAINT_DEFINITION_TAG << " not found,\nDefining equality constraints is optional." << endl;
  }
  else cout << "NSDSXML - loadNSDS : no dynamical systems defined, use of LMGC90 tag." << endl;
}

void NSDSXML::loadNonSmoothDynamicalSystem(NonSmoothDynamicalSystem* nsds)
{
  IN("NSDSXML::loadNSDS( NonSmoothDynamicalSystem* nsds )\n");
  string type;
  string tmp;
  xmlNode* node, *ecDsioNode;
  xmlNode* dsDefinitionNode;
  xmlNode* interactionDefinitionNode, *ecDefinitionNode;
  DSXML* dsxml;
  InteractionXML* interactionXML;
  EqualityConstraintXML *ecXML;
  int number, i;
  char num[32];
  map<int, DSXML*>::iterator it;
  map<int, InteractionXML*>::iterator itinter;
  map<int, EqualityConstraintXML*>::iterator itec;

  if (this->NSDSNode != NULL)
  {
    this->setBVP(nsds->isBVP());

    // at first, we check whether we the tag is LMGC90 tag
    if (SiconosDOMTreeTools::findNodeChild((const xmlNode*)this->NSDSNode, LMGC90_NSDS_TAG) == NULL)
    {
      // creation of the DS_Definition node if necessary
      dsDefinitionNode = SiconosDOMTreeTools::findNodeChild((const xmlNode*)this->NSDSNode, DYNAMICAL_SYSTEM_DEFINITION_TAG);
      if (dsDefinitionNode == NULL)
        dsDefinitionNode = xmlNewChild(this->NSDSNode, NULL, (xmlChar*)DYNAMICAL_SYSTEM_DEFINITION_TAG.c_str(), NULL);

      /*
       * now, creation of the DynamicalSystemXML objects
       */
      for (i = 0; i < nsds->getDSVectorSize(); i++)
      {
        if (nsds->getDynamicalSystem(i)->getDynamicalSystemXML() == NULL)
        {
          type = nsds->getDynamicalSystem(i)->getType();
          number = nsds->getDynamicalSystem(i)->getNumber();
          sprintf(num, "%d", number);
          this->definedDSNumbers.push_back(number);

          // verifies if this Dynamical System has a number which not used
          it = this->DSXMLMap.find(number);
          if (it == this->DSXMLMap.end())
          {
            //node = xmlNewChild( dsDefinitionNode, NULL, (xmlChar*)NSDS_DS.c_str(), NULL );
            if (type == LNLDS)
            {
              node = xmlNewChild(dsDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_NON_LINEARDS_TAG.c_str(), NULL);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              dsxml = new LagrangianNLDSXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getDynamicalSystem(i)->setDynamicalSystemXML(dsxml);

              // creation of the DynamicalSystemXML
              static_cast<LagrangianNLDSXML*>(dsxml)->updateDynamicalSystemXML(node, nsds->getDynamicalSystem(i));

              this->DSXMLMap[number] = dsxml;
            }
            else if (type == LTIDS)
            {
              node = xmlNewChild(dsDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_TIME_INVARIANTDS_TAG.c_str(), NULL);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              dsxml = new LagrangianTIDSXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getDynamicalSystem(i)->setDynamicalSystemXML(dsxml);

              // creation of the DynamicalSystemXML
              static_cast<LagrangianTIDSXML*>(dsxml)->updateDynamicalSystemXML(node, nsds->getDynamicalSystem(i));

              this->DSXMLMap[number] = dsxml;
            }
            else if (type == LSDS)
            {
              node = xmlNewChild(dsDefinitionNode, NULL, (xmlChar*)LINEAR_SYSTEMDS_TAG.c_str(), NULL);
              //xmlNewProp( node, (xmlChar*)NSDS_TYPE.c_str(), (xmlChar*)NSDS_LINEARSYSTEM.c_str() );
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              dsxml = new LinearSystemDSXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getDynamicalSystem(i)->setDynamicalSystemXML(dsxml);

              // creation of the DynamicalSystemXML
              static_cast<LinearSystemDSXML*>(dsxml)->updateDynamicalSystemXML(node, nsds->getDynamicalSystem(i));

              this->DSXMLMap[number] = dsxml;

            }
            else if (type == NLSDS)
            {
              node = xmlNewChild(dsDefinitionNode, NULL, (xmlChar*)NON_LINEAR_SYSTEMDS_TAG.c_str(), NULL);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              dsxml = new DSXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getDynamicalSystem(i)->setDynamicalSystemXML(dsxml);

              // creation of the DynamicalSystemXML
              dsxml->updateDynamicalSystemXML(node, nsds->getDynamicalSystem(i));

              this->DSXMLMap[number] = dsxml;
            }
            else
            {
              XMLException::selfThrow("NSDSXML - loadNSDS error : undefined DS type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
            }
          }
          else
          {
            tmp = num;
            XMLException::selfThrow("NSDSXML - loadNSDS | Error : the Dynamical System number : " + tmp + " already exists!");
          }
        }
      }
    }
    else
    {
      // the LMGC90 tag for DS definition is in the XML file
      //  => specific treatments todo
    }


    // creation of the EqualityConstraint_Defintion if necessary
    if (nsds->getEqualityConstraints().size() > 0)
    {
      ecDefinitionNode = SiconosDOMTreeTools::findNodeChild((const xmlNode*)this->NSDSNode, EQUALITYCONSTRAINT_DEFINITION_TAG);
      if (ecDefinitionNode == NULL)
        ecDefinitionNode = xmlNewChild(this->NSDSNode, NULL, (xmlChar*)EQUALITYCONSTRAINT_DEFINITION_TAG.c_str(), NULL);

      for (i = 0; i < nsds->getEqualityConstraints().size(); i++)
      {
        if (nsds->getEqualityConstraint(i)->getEqualityConstraintXML() == NULL)
        {
          number = nsds->getEqualityConstraint(i)->getNumber();
          sprintf(num, "%d", number);
          this->definedEqualityConstraintNumbers.push_back(number);

          //verifies if the EqualityConstraint has been defined before
          itec = this->equalityConstraintXMLMap.find(number);
          if (itec == this->equalityConstraintXMLMap.end())
          {
            if (nsds->getEqualityConstraint(i)->getType() == LINEAREC)
            {
              node = xmlNewChild(ecDefinitionNode, NULL, (xmlChar*)LINEAR_EC_TAG.c_str(), NULL);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              ecXML = new LinearECXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getEqualityConstraint(i)->setEqualityConstraintXML(ecXML);

              // creation of the DynamicalSystemXML
              static_cast<LinearECXML*>(ecXML)->updateEqualityConstraintXML(node, nsds->getEqualityConstraint(i));

              this->equalityConstraintXMLMap[number] = ecXML;
            }
            else if (nsds->getEqualityConstraint(i)->getType() == LINEARTIEC)
            {
              node = xmlNewChild(ecDefinitionNode, NULL, (xmlChar*)LINEAR_TIME_INVARIANT_EC_TAG.c_str(), NULL);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              ecXML = new LinearTIECXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getEqualityConstraint(i)->setEqualityConstraintXML(ecXML);

              // creation of the DynamicalSystemXML
              static_cast<LinearTIECXML*>(ecXML)->updateEqualityConstraintXML(node, nsds->getEqualityConstraint(i));

              this->equalityConstraintXMLMap[number] = ecXML;
            }
            else if (nsds->getEqualityConstraint(i)->getType() == LAGRANGIANEC)
            {
              node = xmlNewChild(ecDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_EC_TAG.c_str(), NULL);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              ecXML = new LagrangianECXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getEqualityConstraint(i)->setEqualityConstraintXML(ecXML);

              // creation of the DynamicalSystemXML
              static_cast<LagrangianECXML*>(ecXML)->updateEqualityConstraintXML(node, nsds->getEqualityConstraint(i));

              this->equalityConstraintXMLMap[number] = ecXML;
            }
            else if (nsds->getEqualityConstraint(i)->getType() == NLINEAREC)
            {
              node = xmlNewChild(ecDefinitionNode, NULL, (xmlChar*)NON_LINEAR_EC_TAG.c_str(), NULL);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              ecXML = new EqualityConstraintXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getEqualityConstraint(i)->setEqualityConstraintXML(ecXML);

              // creation of the DynamicalSystemXML
              ecXML->updateEqualityConstraintXML(node, nsds->getEqualityConstraint(i));

              this->equalityConstraintXMLMap[number] = ecXML;
            }
            else XMLException::selfThrow("NSDSXML - loadNSDS | Error : the EqualityConstraint type : " + nsds->getEqualityConstraint(i)->getType() + " doesn't exist!");

            /*  end of the save : saving the DynamicalSystem linked to this DSInputOutput */
            ecDsioNode = xmlNewChild(node, NULL, (xmlChar*)DSIO_CONCERNED.c_str(), NULL);
            for (int j = 0; j < nsds->getEqualityConstraint(i)->getDSInputOutputs().size(); j++)
            {
              node = xmlNewChild(ecDsioNode, NULL, (xmlChar*)DSINPUTOUTPUT_TAG.c_str(), NULL);
              number = nsds->getEqualityConstraint(i)->getDSInputOutput(j)->getNumber();
              sprintf(num, "%d", number);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
            }
          }
          else
          {
            tmp = num;
            XMLException::selfThrow("NSDSXML - loadNSDS | Error : the EqualityConstraint number : " + tmp + " already exists!");
          }
        }
      }
    }


    // creation of the Interaction_Defintion if necessary
    if (nsds->getInteractionVectorSize() > 0)
    {
      interactionDefinitionNode = SiconosDOMTreeTools::findNodeChild((const xmlNode*)this->NSDSNode, INTERACTION_DEFINITION_TAG);
      if (interactionDefinitionNode == NULL)
        interactionDefinitionNode = xmlNewChild(this->NSDSNode, NULL, (xmlChar*)INTERACTION_DEFINITION_TAG.c_str(), NULL);

      for (i = 0; i < nsds->getInteractionVectorSize(); i++)
      {
        if (nsds->getInteraction(i)->getInteractionXML() == NULL)
        {
          number = nsds->getInteraction(i)->getNumber();
          sprintf(num, "%d", number);
          this->definedInteractionNumbers.push_back(number);

          // verifies if this Dynamical System has a number which not used
          itinter = this->interactionXMLMap.find(number);
          if (itinter == this->interactionXMLMap.end())
          {
            node = xmlNewChild(interactionDefinitionNode, NULL, (xmlChar*)INTERACTION_TAG.c_str(), NULL);
            xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
            interactionXML = new InteractionXML();

            // linkage between the DynamicalSystem and his DSXML
            nsds->getInteraction(i)->setInteractionXML(interactionXML);

            // creation of the DynamicalSystemXML
            interactionXML->updateInteractionXML(node, nsds->getInteraction(i));

            this->interactionXMLMap[number] = interactionXML;
          }
          else
          {
            tmp = num;
            XMLException::selfThrow("NSDSXML - loadNSDS | Error : the Interaction number : " + tmp + " already exists!");
          }
        }
      }
    }
  }
  else XMLException::selfThrow("NSDSXML - loadNSDS( NSDS* nsds ) Error : no NSDSNode is defined.");
  OUT("NSDSXML::loadNSDS( NonSmoothDynamicalSystem* nsds )\n");
}


void NSDSXML::loadDSXML(xmlNode * rootDSNode)
{
  xmlNode *node;
  int number; //Number of a DS
  map<int, DSXML*>::iterator i;
  string type; //Type of DS
  bool isBVP;

  isBVP = this->isBVP();

  node = SiconosDOMTreeTools::findNodeChild((const xmlNode*)rootDSNode/*, NSDS_DS*/);
  if (node == NULL)
  {
    XMLException::selfThrow("NSDSXML - loadDSXML error : at least one " + DYNAMICAL_SYSTEM_TAG + " must be declared.");
  }

  while (node != NULL)
  {
    DSXML *dsxml;

    number = SiconosDOMTreeTools::getIntegerAttributeValue(node, NUMBER_ATTRIBUTE);

    i = this->DSXMLMap.find(number);
    // we can only add a DS if his number is not already defined
    if (i == this->DSXMLMap.end())
    {
      //type = SiconosDOMTreeTools::getStringAttributeValue(node,NSDS_TYPE);
      /*
       * for XML Schema v1.1 :
       */
      type = (char*)node->name;

      this->definedDSNumbers.push_back(number);

      if (type == LAGRANGIAN_NON_LINEARDS_TAG)
      {
        dsxml = new LagrangianNLDSXML((xmlNode *)node, isBVP);
        this->DSXMLMap[number] = dsxml;
      }
      else if (type == LAGRANGIAN_TIME_INVARIANTDS_TAG)
      {
        dsxml = new LagrangianTIDSXML((xmlNode *)node, isBVP);
        this->DSXMLMap[number] = dsxml;
      }
      else if (type == LINEAR_SYSTEMDS_TAG)
      {
        dsxml = new LinearSystemDSXML((xmlNode *)node, isBVP);
        this->DSXMLMap[number] = dsxml;
      }
      else if (type == NON_LINEAR_SYSTEMDS_TAG)
      {
        dsxml = new DSXML((xmlNode *)node, isBVP);
        this->DSXMLMap[number] = dsxml;
      }
      else
      {
        XMLException::selfThrow("NSDSXML - loadDSXML error : undefined DS type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
      }
    }
    else
    {
      XMLException::selfThrow("NSDSXML - loadDSXML error : wrong DS number : already exists.");
    }

    node = SiconosDOMTreeTools::findFollowNode(node/*, NSDS_DS*/);
  }
}


void NSDSXML::loadInteractionXML(xmlNode * rootInteractionNode)
{
  xmlNode *node;
  int number; //Number of an Interaction
  map<int, InteractionXML*>::iterator i;

  node = SiconosDOMTreeTools::findNodeChild((const xmlNode*)rootInteractionNode, INTERACTION_TAG);

  while (node != NULL)
  {
    InteractionXML *interxml;

    number = SiconosDOMTreeTools::getIntegerAttributeValue(node, NUMBER_ATTRIBUTE);

    cout << "CHARGEMENT INTERACTION number" << number << endl;

    i = interactionXMLMap.find(number);
    if (i == interactionXMLMap.end())
    {
      this->definedInteractionNumbers.push_back(number);
      interxml = new InteractionXML((xmlNode *)node, this->definedDSNumbers);
      this->interactionXMLMap[number] = interxml;
    }
    else
    {
      XMLException::selfThrow("NSDSXML - loadInteractionXML error : wrong INTERACTION number : already exists.");
    }

    node = SiconosDOMTreeTools::findFollowNode(node, INTERACTION_TAG);
  }
}

void NSDSXML::loadEqualityConstraintXML(xmlNode * rootECNode)
{
  xmlNode *node;
  int number; //Number of an EqualityCopnstraint
  map<int, EqualityConstraintXML*>::iterator i;

  node = SiconosDOMTreeTools::findNodeChild((const xmlNode*)rootECNode);

  while (node != NULL)
  {
    EqualityConstraintXML *ecxml;

    number = SiconosDOMTreeTools::getIntegerAttributeValue(node, NUMBER_ATTRIBUTE);

    cout << "CHARGEMENT EqualityConstraint number" << number << endl;

    i = equalityConstraintXMLMap.find(number);
    if (i == equalityConstraintXMLMap.end())
    {
      this->definedEqualityConstraintNumbers.push_back(number);
      ecxml = new EqualityConstraintXML((xmlNode *)node, this->definedDSNumbers);
      this->equalityConstraintXMLMap[number] = ecxml;
    }
    else
    {
      XMLException::selfThrow("NSDSXML - loadEqualityConstraintXML error : wrong EQUALITYCONSTRAINT number : already exists.");
    }

    node = SiconosDOMTreeTools::findFollowNode(node);
  }
}

void NSDSXML::updateNSDSXML(xmlNode* node, NonSmoothDynamicalSystem* nsds)
{
  IN("NSDSXML::updateNSDSXML\n");
  this->NSDSNode = node;
  this->loadNonSmoothDynamicalSystem(nsds);
  OUT("NSDSXML::updateNSDSXML\n");
}

//$Log: NSDSXML.cpp,v $
//Revision 1.47  2005/03/23 15:03:55  jbarbier
//- adaptation to the LMGC90 tags in non smooth dynamical system and strategy
//
//Revision 1.46  2005/03/15 14:44:04  jbarbier
//- pySiconos.i edited to remove local paths
//
//- checkCoherency checks whether the DSInputOutputs and EqualityConstraints have unique numbers
//
//Revision 1.45  2005/03/15 09:57:48  jbarbier
//- EqualityConstraint save OK
//
//Revision 1.44  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.43  2005/03/09 15:30:37  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.42  2005/03/08 12:41:38  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.41  2005/03/04 15:35:27  jbarbier
//- README files added for some samples
//
//- beginning of the refactoring of XML module constants
//
//Revision 1.40  2005/02/04 07:46:21  jbarbier
//- last modification for RollingBalls
//
//Revision 1.39  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.38  2005/01/20 14:44:49  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.37  2004/12/08 12:49:39  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.36  2004/09/16 11:35:25  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.35  2004/09/10 11:26:28  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.34  2004/09/10 08:04:51  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.33  2004/08/23 14:30:03  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.32  2004/08/20 15:26:46  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NSDS and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.31  2004/08/20 07:34:22  jbarbier
//- creation of Model, NSDS in comand program succeed in creating SiconosModelXML,
//NSDSXML
//
//Revision 1.30  2004/08/05 12:44:44  jbarbier
//- loading XML file with no OneStepNSProblem succesfull
//
//- NonLinearSystemDS is now available
//
//Revision 1.29  2004/07/30 14:37:15  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.28  2004/07/29 14:25:44  jbarbier
//- $Log: NSDSXML.cpp,v $
//- Revision 1.47  2005/03/23 15:03:55  jbarbier
//- - adaptation to the LMGC90 tags in non smooth dynamical system and strategy
//-
//- Revision 1.46  2005/03/15 14:44:04  jbarbier
//- - pySiconos.i edited to remove local paths
//-
//- - checkCoherency checks whether the DSInputOutputs and EqualityConstraints have unique numbers
//-
//- Revision 1.45  2005/03/15 09:57:48  jbarbier
//- - EqualityConstraint save OK
//-
//- Revision 1.44  2005/03/14 16:05:27  jbarbier
//- - manual creation of DSInputOutput saving OK
//-
//- - in progress for EqualityConstraint
//-
//- Revision 1.43  2005/03/09 15:30:37  jbarbier
//- - add of LagrangianEC class
//-
//- - in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//-
//- Revision 1.42  2005/03/08 12:41:38  jbarbier
//- - constant variables files modified :
//- Some constants added in SiconosConst
//-
//- all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//-
//- Revision 1.41  2005/03/04 15:35:27  jbarbier
//- - README files added for some samples
//-
//- - beginning of the refactoring of XML module constants
//-
//- Revision 1.40  2005/02/04 07:46:21  jbarbier
//- - last modification for RollingBalls
//-
//- Revision 1.39  2005/01/26 13:50:40  jbarbier
//-
//- - loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//-
//- Revision 1.38  2005/01/20 14:44:49  jbarbier
//- - NSDS class renamed NonSmoothDynamicalSystem
//-
//- - code reduce, some comments remove
//-
//- Revision 1.37  2004/12/08 12:49:39  jbarbier
//- - changes in the XML Schema, respect of the recommandations of the W3C
//- version 1.1
//-
//- - changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//- in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//- for the DS
//-
//- Revision 1.36  2004/09/16 11:35:25  jbarbier
//- - save of the TimeDiscretisation in a XML file in manual creation of the
//- platform which was forgotten is now available.
//-
//- - the save of the platform's data can be done when the platform is created with
//- an XML input file and completed with dynmical systems, interactions, one-step
//- non smooth problem and one-step integrator.
//-
//- Revision 1.35  2004/09/10 11:26:28  charlety
//-
//- _ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//-
//- _ All the tests which worked with the previous version of the vector are OK with the new version.
//-
//- _ Example SICONOS and bouncingBall are OK
//-
//- _ some comments have still to be adapted to NewSiconosVector .
//-
//- _ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//-
//- Revision 1.34  2004/09/10 08:04:51  jbarbier
//- - XML save available for BoundaryCondition and Interaction
//-
//- Revision 1.33  2004/08/23 14:30:03  jbarbier
//- - All the dynamical systems can be created in a comand program and added to a
//- NSDS. The save is OK, but the creation of the boundary conditions is not yet
//- finished.
//-
//- Revision 1.32  2004/08/20 15:26:46  jbarbier
//- - creation of a Model and save in the XML is ok
//- - creation of a NSDS and save in the XML is ok
//- - creation of a NonLinearSystemDS and save in the XML is OK
//-
//- Revision 1.31  2004/08/20 07:34:22  jbarbier
//- - creation of Model, NSDS in comand program succeed in creating SiconosModelXML,
//- NSDSXML
//-
//- Revision 1.30  2004/08/05 12:44:44  jbarbier
//- - loading XML file with no OneStepNSProblem succesfull
//-
//- - NonLinearSystemDS is now available
//-
//- Revision 1.29  2004/07/30 14:37:15  jbarbier
//- - saving methods for DynamicalSystemXML and LagrangianNLDSXML
//- and $Id: NSDSXML.cpp,v 1.47 2005/03/23 15:03:55 jbarbier Exp $ added
//
