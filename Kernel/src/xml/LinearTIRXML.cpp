
#include "LinearTIRXML.h"

LinearTIRXML::LinearTIRXML(): RelationXML()
{
  this->CNode = NULL;
  this->DNode = NULL;
  this->ENode = NULL;
  this->aNode = NULL;
}

LinearTIRXML::LinearTIRXML(xmlNode * LTIRelationNode)
  : RelationXML(LTIRelationNode)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_C)) != NULL)
  {
    this->CNode = node;
  }
  else
  {
    XMLException::selfThrow("LTIRelationXML - constructor error : tag " + LTIR_C + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_D)) != NULL)
  {
    this->DNode = node;
  }
  else
  {
    XMLException::selfThrow("LTIRelationXML - constructor error : tag " + LTIR_D + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_E)) != NULL)
  {
    this->ENode = node;
  }
  else
  {
    XMLException::selfThrow("LTIRelationXML - constructor error : tag " + LTIR_E + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_A)) != NULL)
  {
    this->aNode = node;
  }
  else
  {
    XMLException::selfThrow("LTIRelationXML - constructor error : tag " + LTIR_A + " not found.");
  }
}

LinearTIRXML::~LinearTIRXML()
{
}

void LinearTIRXML::setC(SiconosMatrix *matrix)
{
  if (this->CNode == NULL)
  {
    this->CNode = SiconosDOMTreeTools::createMatrixNode(this->rootRelationXMLNode, LTIR_C, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixValue(CNode, matrix);
}

void LinearTIRXML::setD(SiconosMatrix *matrix)
{
  if (this->DNode == NULL)
  {
    this->DNode = SiconosDOMTreeTools::createMatrixNode(this->rootRelationXMLNode, LTIR_D, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixValue(DNode, matrix);
}

void LinearTIRXML::setE(SiconosMatrix *matrix)
{
  if (this->ENode == NULL)
  {
    this->ENode = SiconosDOMTreeTools::createMatrixNode(this->rootRelationXMLNode, LTIR_E, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixValue(ENode, matrix);
}

void LinearTIRXML::setA(SiconosVector *vector)
{
  if (this->aNode == NULL)
  {
    this->aNode = SiconosDOMTreeTools::createVectorNode(this->rootRelationXMLNode, LTIR_A, vector);
  }
  else SiconosDOMTreeTools::setSiconosVectorValue(aNode, vector);
}

