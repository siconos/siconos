#include "OneStepIntegratorXML.h"
using namespace std;

OneStepIntegratorXML::OneStepIntegratorXML():
  rootIntegratorXMLNode(NULL), rNode(NULL), DSConcernedNode(NULL)
{}

OneStepIntegratorXML::OneStepIntegratorXML(xmlNode * OneStepIntegratorNode, map<int, bool> definedDSNumbers):
  rootIntegratorXMLNode(OneStepIntegratorNode), rNode(NULL), DSConcernedNode(NULL)
{
  xmlNode * node;
  if ((node = SiconosDOMTreeTools::findNodeChild(OneStepIntegratorNode, OSI_DS_CONCERNED)) != NULL)
  {
    DSConcernedNode = node;
    loadOneStepIntegratorConcernedDS(node, definedDSNumbers);
  }
  else
  {
    XMLException::selfThrow("OneStepIntegratorXML - Constructor error : tag " + OSI_DS_CONCERNED + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(OneStepIntegratorNode, OSI_R)) != NULL)
    rNode = node;
  else
  {
    //XMLException::selfThrow("OneStepIntegratorXML - Constructor error : tag " + OSI_R + " not found.");
    cout << "OneStepIntegratorXML - Constructor : Warning : tag " << OSI_R << " not found." << endl;
    rNode = NULL;
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

