
#include "OneStepIntegratorXML.h"

#include "check.h"


OneStepIntegratorXML::OneStepIntegratorXML()
{
  this->rootIntegratorXMLNode = NULL;
  this->rNode = NULL;
  this->DSConcernedNode = NULL;
}

OneStepIntegratorXML::OneStepIntegratorXML(xmlNode * OneStepIntegratorNode, map<int, bool> definedDSNumbers)
{
  xmlNode * node;
  this->rootIntegratorXMLNode = OneStepIntegratorNode;
  if ((node = SiconosDOMTreeTools::findNodeChild(OneStepIntegratorNode, OSI_DS_CONCERNED)) != NULL)
  {
    this->DSConcernedNode = node;
    this->loadOneStepIntegratorConcernedDS(node, definedDSNumbers);
  }
  else
  {
    XMLException::selfThrow("OneStepIntegratorXML - Constructor error : tag " + OSI_DS_CONCERNED + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(OneStepIntegratorNode, OSI_R)) != NULL)
  {
    this->rNode = node;
  }
  else
  {
    //XMLException::selfThrow("OneStepIntegratorXML - Constructor error : tag " + OSI_R + " not found.");
    cout << "OneStepIntegratorXML - Constructor : Warning : tag " << OSI_R << " not found." << endl;
    this->rNode = NULL;
  }
}

OneStepIntegratorXML::~OneStepIntegratorXML()
{
}


void OneStepIntegratorXML::loadOneStepIntegratorConcernedDS(xmlNode * DSConcernedNode, map<int, bool> definedDSNumbers)
{
  xmlNode *DSnode;
  int number;
  //int size = SiconosDOMTreeTools::getIntegerAttributeValue(DSConcernedNode, OSI_SIZE);
  int size = 0;
  int i = 0;
  map<int, bool>::iterator iter;

  DSNumbersVector.clear();

  if ((DSnode = SiconosDOMTreeTools::findNodeChild((const xmlNode*)DSConcernedNode, DYNAMICAL_SYSTEM_TAG)) == NULL)
    XMLException::selfThrow("OneStepIntegratorXML - loadOneStepIntegratonConcernedDS error : at least one " + DYNAMICAL_SYSTEM_TAG + " tag must be declared in " + OSI_DS_CONCERNED + " tag.");

  size = SiconosDOMTreeTools::getNodeChildrenNumber(DSConcernedNode);
  //  cout<<"number of children = "<<size<<endl;
  //  cout<<"# Press <<ENTER>> !!"<<endl;
  //  getchar();

  while ((DSnode != NULL) && (i < size))
  {
    number = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, NUMBER_ATTRIBUTE);
    //cout<<"NUMBER "<<number<<endl; getchar();
    //Verifying DS number exists and not used by another OneStepIntegrator
    iter = definedDSNumbers.find(number);
    if (iter == definedDSNumbers.end())
    {
      char errorMsg[1024];
      sprintf(errorMsg, "OneStepIntegratorXML - loadOneStepIntegratorConcernedDS error : in a tag %s you define a DS number who doesn't exist : %d.", OSI_DS_CONCERNED.c_str(), number);
      XMLException::selfThrow(errorMsg);

    }

    if (definedDSNumbers[number] == false)
    {
      char errorMsg[1024];
      sprintf(errorMsg, "OneStepIntegratorXML - loadOneStepIntegratorConcernedDS error : in tag %s you define a DS number : %d who is already used in another OneStepIntegrator tag.", OSI_DS_CONCERNED.c_str(), number);
      XMLException::selfThrow(errorMsg);
    }

    // Now DS number is used in a OneStepIntegrator definition
    definedDSNumbers[number] = true;

    //this->DSNumbersVector[i]=number;
    this->DSNumbersVector.push_back(number);

    DSnode = SiconosDOMTreeTools::findFollowNode(DSnode, DYNAMICAL_SYSTEM_TAG);

    i++;
  }
  if (i < size)
  {
    XMLException::selfThrow("OneStepIntegratorXML - loadOneStepIntegratorConcernedDS error : the size attribute given in the tag " + OSI_DS_CONCERNED + " is not correct.");
  }
}


void OneStepIntegratorXML::updateOneStepIntegratorXML(xmlNode* node, OneStepIntegrator* osi)
{
  IN("OneStepIntegratorXML::updateOneStepIntegratorXML\n");
  this->rootIntegratorXMLNode = node;
  //  this->loadOneStepIntegrator( osi );
  OUT("OneStepIntegratorXML::updateOneStepIntegratorXML\n");
}

void OneStepIntegratorXML::setDSConcerned(vector<int>* ds)
{
  if (this->DSConcernedNode == NULL)
  {
    this->DSNumbersVector = *ds;

    //save in the DOM tree
    char num[32];
    xmlNode* node;
    /*
     * creation of the DS_Concerned node
     */
    this->DSConcernedNode = xmlNewChild(this->rootIntegratorXMLNode, NULL, (xmlChar*)OSI_DS_CONCERNED.c_str(), NULL);
    sprintf(num, "%i", this->DSNumbersVector.size());
    xmlNewProp(this->DSConcernedNode, (xmlChar*)SIZE_ATTRIBUTE.c_str(), (xmlChar*)num);

    for (int i = 0; i < this->DSNumbersVector.size(); i++)
    {
      node = xmlNewChild(this->DSConcernedNode, NULL, (xmlChar*)DYNAMICAL_SYSTEM_TAG.c_str(), NULL);
      sprintf(num, "%i", this->DSNumbersVector[i]);
      xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
    }
  }
  else
  {
    /* \todo : when the user complete the data of the XML input file
     * it must check already defined dynamical system number before adding them */
  }
}

//$Log: OneStepIntegratorXML.cpp,v $
//Revision 1.24  2005/03/08 14:23:45  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.23  2005/01/10 17:06:37  jbarbier
//- attribute "size" is now unused in the code
//
//- xml schema v1.2 is in progress
//
//Revision 1.22  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.21  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.20  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.19  2004/09/14 13:49:58  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.18  2004/09/10 11:26:28  charlety
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
//Revision 1.17  2004/07/29 14:25:44  jbarbier
