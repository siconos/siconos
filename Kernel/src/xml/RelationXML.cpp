
#include "RelationXML.h"


RelationXML::RelationXML()
{
  this->computeInputNode = NULL;
  this->computeOutputNode = NULL;
}

RelationXML::RelationXML(xmlNode *relationNode)
{
  IN("RelationXML::RelationXML(xmlNode*)\n");
  xmlNode *node;
  string type((char*)relationNode->name);
  this->rootRelationXMLNode = relationNode;

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootRelationXMLNode, COMPUTE_INPUT_TAG)) != NULL)
  {
    this->computeInputNode = node;
  }
  else
  {
    //if( SiconosDOMTreeTools::getStringAttributeValue(this->rootRelationXMLNode, RELATION_TYPE) == RELATION_LNL )
    if (LAGRANGIAN_NON_LINEAR_RELATION_TAG == type)
    {
      XMLException::selfThrow("RelationXML - RelationXML::RelationXML(xmlNode *relationNode) error : tag " + COMPUTE_INPUT_TAG + " not found.");
    }
    //else cout<< "RelationXML - RelationXML::RelationXML(xmlNode *relationNode) Warning : tag " << RELATION_INPUT << " not found. This is attribute is optional for this relation" <<endl;
    this->computeInputNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootRelationXMLNode, COMPUTE_OUTPUT_TAG)) != NULL)
  {
    this->computeOutputNode = node;
  }
  else
  {
    //if( SiconosDOMTreeTools::getStringAttributeValue(this->rootRelationXMLNode, RELATION_TYPE) == RELATION_LNL )
    if (LAGRANGIAN_NON_LINEAR_RELATION_TAG == type)
    {
      XMLException::selfThrow("RelationXML - RelationXML::RelationXML(xmlNode *relationNode) error : tag " + COMPUTE_OUTPUT_TAG + " not found.");
    }
    //else cout<< "RelationXML - RelationXML::RelationXML(xmlNode *relationNode) Warning : tag " << RELATION_OUTPUT << " not found. This is attribute is optional for this relation" <<endl;
    this->computeOutputNode = NULL;
  }
  OUT("RelationXML::RelationXML(xmlNode*)\n");
}

RelationXML::~RelationXML()
{}

void RelationXML::updateRelationXML(xmlNode* node, Relation* rel)
{
  IN("RelationXML::updateRelationXML\n");
  this->rootRelationXMLNode = node;
  OUT("RelationXML::updateRelationXML\n");
}

