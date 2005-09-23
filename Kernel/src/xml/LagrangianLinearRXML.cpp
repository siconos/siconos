#include "LagrangianLinearRXML.h"
using namespace std;

LagrangianLinearRXML::LagrangianLinearRXML():
  LagrangianRXML(), HNode(NULL), bNode(NULL)
{}

LagrangianLinearRXML::LagrangianLinearRXML(xmlNode * LLRelationNode)
  : LagrangianRXML(LLRelationNode), HNode(NULL), bNode(NULL)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, LLR_H)) != NULL)
    HNode = node;
  else
    XMLException::selfThrow("LLRelationXML - constructor error : tag " + LLR_H + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, LLR_B)) != NULL)
    bNode = node;
  else
    XMLException::selfThrow("LLRelationXML - constructor error : tag " + LLR_B + " not found.");
}

LagrangianLinearRXML::~LagrangianLinearRXML()
{}

void LagrangianLinearRXML::setH(const SiconosMatrix& matrix)
{
  if (HNode == NULL)
  {
    HNode = SiconosDOMTreeTools::createMatrixNode(rootRelationXMLNode, LLR_H, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(HNode, matrix);
}

void LagrangianLinearRXML::setB(const SiconosVector &vec)
{
  if (bNode == NULL)
  {
    bNode = SiconosDOMTreeTools::createVectorNode(rootRelationXMLNode, LLR_B, vec);
  }
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(bNode, vec);
}
