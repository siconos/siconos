
#include "NSDSXML.h"

#include "LagrangianDSXML.h"
#include "LagrangianLinearTIDSXML.h"
#include "LinearSystemDSXML.h"

#include "LinearECXML.h"
#include "LinearTIECXML.h"
#include "LagrangianECXML.h"
#include "LagrangianLinearECXML.h"
#include "LinearDSIOXML.h"
#include "LagrangianDSIOXML.h"
#include "LagrangianLinearDSIOXML.h"

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
  xmlNode *childNode;

  if ((node = SiconosDOMTreeTools::findNodeChild(NSDSNode, LMGC90_NSDS_TAG)) == NULL)
  {
    // at first, we load the DSInputOutputs because we need them to load properly the DSXML
    if ((node = SiconosDOMTreeTools::findNodeChild(NSDSNode, DSINPUTOUTPUT_DEFINITION_TAG)) != NULL)
      this->loadDSInputOutputXML(node);
    else
      cout << "DSXML - loadNSDS WARNING : tag " << DSINPUTOUTPUT_DEFINITION_TAG << " not found,\nDefining DS InputOutput is optional." << endl;

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
  IN("NSDSXML::loadNonSmoothDynamicalSystem( NonSmoothDynamicalSystem* nsds )\n");
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
        if (nsds->getDynamicalSystem(i)->getDynamicalSystemXMLPtr() == NULL)
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
              dsxml = new LagrangianDSXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getDynamicalSystem(i)->setDynamicalSystemXMLPtr(dsxml);

              // creation of the DynamicalSystemXML
              static_cast<LagrangianDSXML*>(dsxml)->updateDynamicalSystemXML(node, nsds->getDynamicalSystem(i));

              this->DSXMLMap[number] = dsxml;
            }
            else if (type == LTIDS)
            {
              node = xmlNewChild(dsDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_TIME_INVARIANTDS_TAG.c_str(), NULL);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              dsxml = new LagrangianLinearTIDSXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getDynamicalSystem(i)->setDynamicalSystemXMLPtr(dsxml);

              // creation of the DynamicalSystemXML
              static_cast<LagrangianLinearTIDSXML*>(dsxml)->updateDynamicalSystemXML(node, nsds->getDynamicalSystem(i));

              this->DSXMLMap[number] = dsxml;
            }
            else if (type == LDS)
            {
              node = xmlNewChild(dsDefinitionNode, NULL, (xmlChar*)LINEAR_SYSTEMDS_TAG.c_str(), NULL);
              //xmlNewProp( node, (xmlChar*)NSDS_TYPE.c_str(), (xmlChar*)NSDS_LINEARSYSTEM.c_str() );
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              dsxml = new LinearSystemDSXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getDynamicalSystem(i)->setDynamicalSystemXMLPtr(dsxml);

              // creation of the DynamicalSystemXML
              static_cast<LinearSystemDSXML*>(dsxml)->updateDynamicalSystemXML(node, nsds->getDynamicalSystem(i));

              this->DSXMLMap[number] = dsxml;

            }
            else if (type == NLDS)
            {
              node = xmlNewChild(dsDefinitionNode, NULL, (xmlChar*)NON_LINEAR_SYSTEMDS_TAG.c_str(), NULL);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              dsxml = new DSXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getDynamicalSystem(i)->setDynamicalSystemXMLPtr(dsxml);

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
        else
          cout << "## /!\\ the DynamicalSystem : " << nsds->getDynamicalSystem(i)->getType() << " number " << nsds->getDynamicalSystem(i)->getNumber() <<
               ", has already an XML object." << endl;
      }
    }
    else
    {
      // the LMGC90 tag for DS definition is in the XML file
      // \todo !!!!!!  => specific treatments todo
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
            else if (nsds->getEqualityConstraint(i)->getType() == LAGRANGIANLINEAREC)
            {
              node = xmlNewChild(ecDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_LINEAR_EC_TAG.c_str(), NULL);
              xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
              ecXML = new LagrangianECXML();

              // linkage between the DynamicalSystem and his DSXML
              nsds->getEqualityConstraint(i)->setEqualityConstraintXML(ecXML);

              // creation of the DynamicalSystemXML
              static_cast<LagrangianLinearECXML*>(ecXML)->updateEqualityConstraintXML(node, nsds->getEqualityConstraint(i));

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
  else XMLException::selfThrow("NSDSXML - loadNonSmoothDynamicalSystem( NSDS* nsds ) Error : no NSDSNode is defined.");
  OUT("NSDSXML::loadNonSmoothDynamicalSystem( NonSmoothDynamicalSystem* nsds )\n");
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
      type = (char*)node->name;
      cout << "#~#~ Loading DynamicalSystem number : " << number << " - " << type << endl;
      this->definedDSNumbers.push_back(number);

      if (type == LAGRANGIAN_NON_LINEARDS_TAG)
      {
        dsxml = new LagrangianDSXML((xmlNode *)node, isBVP);
        this->DSXMLMap[number] = dsxml;
        dsxml->setDSInputOutputXML(this->getDSInputOutputXMLRelatingToDS(number));
      }
      else if (type == LAGRANGIAN_TIME_INVARIANTDS_TAG)
      {
        dsxml = new LagrangianLinearTIDSXML((xmlNode *)node, isBVP);
        this->DSXMLMap[number] = dsxml;
        dsxml->setDSInputOutputXML(this->getDSInputOutputXMLRelatingToDS(number));
      }
      else if (type == LINEAR_SYSTEMDS_TAG)
      {
        dsxml = new LinearSystemDSXML((xmlNode *)node, isBVP);
        this->DSXMLMap[number] = dsxml;
        dsxml->setDSInputOutputXML(this->getDSInputOutputXMLRelatingToDS(number));
      }
      else if (type == NON_LINEAR_SYSTEMDS_TAG)
      {
        dsxml = new DSXML((xmlNode *)node, isBVP);
        this->DSXMLMap[number] = dsxml;
        dsxml->setDSInputOutputXML(this->getDSInputOutputXMLRelatingToDS(number));
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

map<int, DSInputOutputXML*> NSDSXML::getDSInputOutputXMLRelatingToDS(int number)
{
  int i;
  map<int, DSInputOutputXML*> m;
  vector<int> v;

  map<int, DSInputOutputXML*>::iterator iter;
  for (iter = dsInputOutputXMLMap.begin(); iter != dsInputOutputXMLMap.end(); iter++)
  {
    v = (*iter).second->getDSConcerned();
    for (i = 0; i < v.size(); i++)
    {
      //      cout<<"** NSDSXML::getDSInputOutputXMLRelatingToDS v["<<i<<"] == "<<v[i]<<endl;
      if (v[i] == number)
      {
        m[(*iter).first] = (*iter).second;
        //        cout<<"** NSDSXML::getDSInputOutputXMLRelatingToDS ==> "<<(*iter).first<<" - "<<(*iter).second->getType()<<endl;
      }
    }
  }

  return m;
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

void NSDSXML::loadDSInputOutputXML(xmlNode * rootdsioNode)
{
  xmlNode *node;
  int number;
  string type;
  map<int, DSInputOutputXML*>::iterator i;
  DSInputOutputXML *dsioxml;

  node = SiconosDOMTreeTools::findNodeChild((const xmlNode*)rootdsioNode);


  while (node != NULL)
  {
    number = SiconosDOMTreeTools::getIntegerAttributeValue(node, NUMBER_ATTRIBUTE);
    type = (char*)node->name;
    cout << "#~#~ Loading DSInputOutput number : " << number << " - " << type << endl;

    i = dsInputOutputXMLMap.find(number);
    if (i == dsInputOutputXMLMap.end())
    {
      if (type == LINEAR_DSIO_TAG)
      {
        this->definedDSInputOutputNumbers.push_back(number);
        dsioxml = new LinearDSIOXML((xmlNode *)node/*, this->definedDSNumbers*/);
        this->dsInputOutputXMLMap[number] = dsioxml;
      }
      else if (type == NON_LINEAR_DSIO_TAG)
      {
        this->definedDSInputOutputNumbers.push_back(number);
        dsioxml = new DSInputOutputXML((xmlNode *)node/*, this->definedDSNumbers*/);
        this->dsInputOutputXMLMap[number] = dsioxml;
      }
      else if (type == LAGRANGIAN_DSIO_TAG)
      {
        this->definedDSInputOutputNumbers.push_back(number);
        dsioxml = new LagrangianDSIOXML((xmlNode *)node/*, this->definedDSNumbers*/);
        this->dsInputOutputXMLMap[number] = dsioxml;
      }
      else if (type == LAGRANGIAN_LINEAR_DSIO_TAG)
      {
        this->definedDSInputOutputNumbers.push_back(number);
        dsioxml = new LagrangianLinearDSIOXML((xmlNode *)node/*, this->definedDSNumbers*/);
        this->dsInputOutputXMLMap[number] = dsioxml;
      }
      else
        XMLException::selfThrow("NSDSXML - loadDSInputOutputXML error : wrong DSInputOutput number : already exists.");
    }
    else
    {
      XMLException::selfThrow("NSDSXML - loadDSInputOutputXML error : wrong DSINPUTOUTPUT number : already exists.");
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

