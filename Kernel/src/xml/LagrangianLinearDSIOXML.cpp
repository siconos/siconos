#include "LagrangianLinearDSIOXML.h"
using namespace std;

LagrangianLinearDSIOXML::LagrangianLinearDSIOXML(): LagrangianDSIOXML()
{
  this->HNode = NULL;
  this->bNode = NULL;
}

LagrangianLinearDSIOXML::LagrangianLinearDSIOXML(xmlNode * dsioNode)
// : DSInputOutputXML(dsioNode)
{
  xmlNode *node;
  this->rootDSIOXMLNode = dsioNode;
  if ((node = SiconosDOMTreeTools::findNodeChild(dsioNode, LLDSIO_H)) != NULL)
    this->HNode = node;
  else
    XMLException::selfThrow("LagrangianLinearDSIOXML - constructor error : tag " + LLDSIO_H + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(dsioNode, LLDSIO_B)) != NULL)
    this->bNode = node;
  else
    XMLException::selfThrow("LagrangianLinearDSIOXML - constructor error : tag " + LLDSIO_B + " not found.");
}

LagrangianLinearDSIOXML::~LagrangianLinearDSIOXML()
{}

void LagrangianLinearDSIOXML::setH(SiconosMatrix *matrix)
{
  if (this->HNode == NULL)
  {
    this->HNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSIOXMLNode, LLDSIO_H, *matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(this->HNode, *matrix);
}

void LagrangianLinearDSIOXML::setB(SiconosVector *vector)
{
  if (this->bNode == NULL)
  {
    this->bNode = SiconosDOMTreeTools::createVectorNode(this->rootDSIOXMLNode, LLDSIO_B, *vector);
  }
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->bNode, *vector);
}
