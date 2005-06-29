#include "RelationXML.h"
using namespace std;

RelationXML::RelationXML():
  rootRelationXMLNode(NULL), computeInputNode(NULL), computeOutputNode(NULL)
{}

RelationXML::RelationXML(xmlNode *relationNode):
  rootRelationXMLNode(relationNode), computeInputNode(NULL), computeOutputNode(NULL)
{
  IN("RelationXML::RelationXML(xmlNode*)\n");
  xmlNode *node;
  string type((char*)relationNode->name);
  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootRelationXMLNode, COMPUTE_INPUT_TAG)) != NULL)
    computeInputNode = node;
  else
  {
    if (LAGRANGIAN_NON_LINEAR_RELATION_TAG == type)
      XMLException::selfThrow("RelationXML - RelationXML::RelationXML(xmlNode *relationNode) error : tag " + COMPUTE_INPUT_TAG + " not found.");
  }
  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootRelationXMLNode, COMPUTE_OUTPUT_TAG)) != NULL)
    computeOutputNode = node;
  else
  {
    if (LAGRANGIAN_NON_LINEAR_RELATION_TAG == type)
      XMLException::selfThrow("RelationXML - RelationXML::RelationXML(xmlNode *relationNode) error : tag " + COMPUTE_OUTPUT_TAG + " not found.");
  }
  OUT("RelationXML::RelationXML(xmlNode*)\n");
}

RelationXML::~RelationXML()
{}


void RelationXML::setComputeInputPlugin(const string&  plugin)
{
  if (computeInputNode == NULL)
  {
    computeInputNode = SiconosDOMTreeTools::createSingleNode(rootRelationXMLNode, COMPUTE_INPUT_TAG);
    xmlNewProp(computeInputNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
  }
  else SiconosDOMTreeTools::setStringAttributeValue(computeInputNode, PLUGIN_ATTRIBUTE, plugin);
}

void RelationXML::setComputeOutputPlugin(const string&  plugin)
{
  if (computeOutputNode == NULL)
  {
    computeOutputNode = SiconosDOMTreeTools::createSingleNode(rootRelationXMLNode, COMPUTE_OUTPUT_TAG);
    xmlNewProp(computeOutputNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
  }
  else SiconosDOMTreeTools::setStringAttributeValue(computeOutputNode, PLUGIN_ATTRIBUTE, plugin);
}

void RelationXML::updateRelationXML(xmlNode* node, Relation* rel)
{
  IN("RelationXML::updateRelationXML\n");
  rootRelationXMLNode = node;
  OUT("RelationXML::updateRelationXML\n");
}

string RelationXML::getComputeInputPlugin() const
{
  if (!isComputeInputPlugin())
    XMLException::selfThrow("RelationXML - getComputeInputPlugin : computeInput is not calculated from a plugin");
  return  SiconosDOMTreeTools::getStringAttributeValue(computeInputNode, PLUGIN_ATTRIBUTE);
}
string RelationXML::getComputeOutputPlugin() const
{
  if (!isComputeOutputPlugin())
    XMLException::selfThrow("RelationXML - getComputeOutputPlugin : computeOutput is not calculated from a plugin");
  return  SiconosDOMTreeTools::getStringAttributeValue(computeOutputNode, PLUGIN_ATTRIBUTE);
}
