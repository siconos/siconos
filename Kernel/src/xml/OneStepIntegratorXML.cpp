#include "OneStepIntegratorXML.h"
using namespace std;

OneStepIntegratorXML::OneStepIntegratorXML():
  rNode(NULL), rootIntegratorXMLNode(NULL), DSConcernedNode(NULL)
{}

OneStepIntegratorXML::OneStepIntegratorXML(xmlNode * OneStepIntegratorNode, map<int, bool> definedDSNumbers):
  rNode(NULL), rootIntegratorXMLNode(OneStepIntegratorNode), DSConcernedNode(NULL)
{
  xmlNode * node;
  if ((node = SiconosDOMTreeTools::findNodeChild(OneStepIntegratorNode, OSI_DS_CONCERNED)) != NULL)
  {
    DSConcernedNode = node;
    loadOneStepIntegratorConcernedDS(node, definedDSNumbers);
  }
  else
    XMLException::selfThrow("OneStepIntegratorXML - Constructor error : tag " + OSI_DS_CONCERNED + " not found.");
  if ((node = SiconosDOMTreeTools::findNodeChild(OneStepIntegratorNode, OSI_R)) != NULL)
    rNode = node;
}

OneStepIntegratorXML::~OneStepIntegratorXML()
{}

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
  while ((DSnode != NULL) && (i < size))
  {
    number = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, NUMBER_ATTRIBUTE);
    //Verifying DS number exists and not used by another OneStepIntegrator
    iter = definedDSNumbers.find(number);
    if (iter == definedDSNumbers.end())
      XMLException::selfThrow("OneStepIntegratorXML - loadOneStepIntegratorConcernedDS error : in tag" + OSI_DS_CONCERNED + "undefined DS number");

    if (definedDSNumbers[number] == false)
      XMLException::selfThrow("OneStepIntegratorXML - loadOneStepIntegratorConcernedDS error : in tag" + OSI_DS_CONCERNED + "already used DS number");

    // Now DS number is used in a OneStepIntegrator definition
    definedDSNumbers[number] = true;

    DSNumbersVector.push_back(number);

    DSnode = SiconosDOMTreeTools::findFollowNode(DSnode, DYNAMICAL_SYSTEM_TAG);
    i++;
  }
  if (i < size)
    XMLException::selfThrow("OneStepIntegratorXML - loadOneStepIntegratorConcernedDS error : wrong size attribute in the tag " + OSI_DS_CONCERNED);
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
  if (DSConcernedNode == NULL)
  {
    DSNumbersVector = *ds;

    //save in the DOM tree
    char num[32];
    xmlNode* node;
    /*
     * creation of the DS_Concerned node
     */
    DSConcernedNode = xmlNewChild(rootIntegratorXMLNode, NULL, (xmlChar*)OSI_DS_CONCERNED.c_str(), NULL);
    sprintf(num, "%i", DSNumbersVector.size());
    xmlNewProp(DSConcernedNode, (xmlChar*)SIZE_ATTRIBUTE.c_str(), (xmlChar*)num);

    for (unsigned int i = 0; i < DSNumbersVector.size(); i++)
    {
      node = xmlNewChild(DSConcernedNode, NULL, (xmlChar*)DYNAMICAL_SYSTEM_TAG.c_str(), NULL);
      sprintf(num, "%i", DSNumbersVector[i]);
      xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
    }
  }
  else
  {
    /* \todo : when the user complete the data of the XML input file
     * it must check already defined dynamical system number before adding them */
  }
}

