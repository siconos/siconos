
#include "DSInputOutputXML.h"


DSInputOutputXML::DSInputOutputXML()
{
  //  this->computeInputNode = NULL;
  //  this->computeOutputNode = NULL;
  this->HNode = NULL;
}

DSInputOutputXML::DSInputOutputXML(xmlNode *dsioNode/*, vector<int> definedDSNumbers*/)
{
  IN("DSInputOutputXML::DSInputOutputXML(xmlNode*)\n");
  xmlNode *node;
  //string type ( (char*)dsioNode->name );
  this->rootDSIOXMLNode = dsioNode;

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSIOXMLNode, DSINPUTOUTPUT_H)) != NULL)
  {
    this->HNode = node;
  }
  else
  {
    XMLException::selfThrow("DSInputOutputXML - DSInputOutputXML(xmlNode *dsioNode) error : tag " + DSINPUTOUTPUT_H + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSIOXMLNode, DS_CONCERNED)) != NULL)
  {
    this->dsConcernedNode = node;
    loadDSIOConcernedDS(node/*, definedDSNumbers*/);
  }
  else
  {
    XMLException::selfThrow("DSInputOutputXML - DSInputOutputXML(xmlNode *dsioNode) error : tag " + DS_CONCERNED + " not found.");
  }
  OUT("DSInputOutputXML::DSInputOutputXML(xmlNode*)\n");
}

DSInputOutputXML::~DSInputOutputXML()
{}


void DSInputOutputXML::loadDSIOConcernedDS(xmlNode * DSConcernedNode/*, vector<int> definedDSNumbers*/)
{
  IN("DSInputOutputXML::loadDSIOConcernedDS\n ");

  xmlNode *DSnode;
  int number;
  //int size = SiconosDOMTreeTools::getIntegerAttributeValue(DSConcernedNode, INTERACTION_SIZE);
  int size = 0;
  int i = 0, j = 0;

  if ((DSnode = SiconosDOMTreeTools::findNodeChild((const xmlNode*)DSConcernedNode, DYNAMICAL_SYSTEM_TAG)) == NULL)
  {
    XMLException::selfThrow("DSInputOutputXML - loadDSIOConcernedDS error : at least one couple of " + DYNAMICAL_SYSTEM_TAG + " must be declared in " + DS_CONCERNED + " tag.");
  }

  size = SiconosDOMTreeTools::getNodeChildrenNumber(DSConcernedNode);
  while ((DSnode != NULL) && (i < size))
  {
    number = SiconosDOMTreeTools::getIntegerAttributeValue(DSnode, NUMBER_ATTRIBUTE);

    // \todo : verifying that the DS are defined before
    //    j = 0;
    //    //Verifying that the DS number exists
    //    while ((j<definedDSNumbers.size()) && (definedDSNumbers[j]!=number)) {  j++; }
    //
    //    if (j==definedDSNumbers.size())
    //    {
    //      char errorMsg[1024];
    //      sprintf(errorMsg, "DSInputOutputXML - loadDSInputOutputConcernedDS error : in a tag %s you define couple of DS with a DS number who doesn't exist : %d.", INTERACTION_DS_CONCERNED.c_str(), number1);
    //      XMLException::selfThrow(errorMsg);
    //    }

    this->definedDSNumbers.push_back(number);

    DSnode = SiconosDOMTreeTools::findFollowNode(DSnode, DYNAMICAL_SYSTEM_TAG);

    i++;
  }

  OUT("DSInputOutputXML::loadDSIOConcernedDS\n ");
}

void DSInputOutputXML::updateDSInputOutputXML(xmlNode* node, DSInputOutput* dsio)
{
  IN("DSInputOutputXML::updateDSInputOutputXML\n");
  this->rootDSIOXMLNode = node;
  OUT("DSInputOutputXML::updateDSInputOutputXML\n");
}

void DSInputOutputXML::setDSConcerned(vector<int> dsConcerned)
{
  if (this->dsConcernedNode == NULL)
  {
    int i;
    int number;
    xmlNode* node;

    //save in the DOM tree
    char num[32];

    /*
     * creation of the DS_Concerned node
     */
    this->dsConcernedNode = xmlNewChild(this->rootDSIOXMLNode, NULL, (xmlChar*)DS_CONCERNED.c_str(), NULL);

    for (i = 0; i < this->definedDSNumbers.size(); i++)
    {
      node = xmlNewChild(this->dsConcernedNode, NULL, (xmlChar*)DYNAMICAL_SYSTEM_TAG.c_str(), NULL);
      sprintf(num, "%i", this->definedDSNumbers[i]);
      xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
    }
  }
  else
  {
    /* \todo : when DSIO have been given in the XML input file and that the user has added new ones */
  }
}

