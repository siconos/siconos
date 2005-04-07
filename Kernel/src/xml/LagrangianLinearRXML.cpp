
#include "LagrangianLinearRXML.h"

LagrangianLinearRXML::LagrangianLinearRXML(): RelationXML()
{
  this->HNode = NULL;
  this->bNode = NULL;
}

LagrangianLinearRXML::LagrangianLinearRXML(xmlNode * LLRelationNode)
  : RelationXML(LLRelationNode)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, LLR_H)) != NULL)
  {
    this->HNode = node;
  }
  else
  {
    XMLException::selfThrow("LLRelationXML - constructor error : tag " + LLR_H + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, LLR_B)) != NULL)
  {
    this->bNode = node;
  }
  else
  {
    XMLException::selfThrow("LLRelationXML - constructor error : tag " + LLR_B + " not found.");
  }
}

LagrangianLinearRXML::~LagrangianLinearRXML()
{}

void LagrangianLinearRXML::setH(SiconosMatrix *matrix)
{
  if (this->HNode == NULL)
  {
    this->HNode = SiconosDOMTreeTools::createMatrixNode(this->rootRelationXMLNode, LLR_H, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixValue(this->HNode, matrix);
}

void LagrangianLinearRXML::setB(SiconosVector *vector)
{
  if (this->bNode == NULL)
  {
    this->bNode = SiconosDOMTreeTools::createVectorNode(this->rootRelationXMLNode, LLR_B, vector);
  }
  else SiconosDOMTreeTools::setSiconosVectorValue(this->bNode, vector);
}
