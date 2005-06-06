#include "LinearTIRXML.h"
using namespace std;

LinearTIRXML::LinearTIRXML():
  RelationXML(), CNode(NULL), DNode(NULL), ENode(NULL), aNode(NULL)
{}

LinearTIRXML::LinearTIRXML(xmlNode * LTIRelationNode)
  : RelationXML(LTIRelationNode), CNode(NULL), DNode(NULL), ENode(NULL), aNode(NULL)
{
  xmlNode *node;
  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_C)) != NULL)
    CNode = node;
  else
    XMLException::selfThrow("LTIRelationXML - constructor error : tag " + LTIR_C + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_D)) != NULL)
    DNode = node;
  else
    XMLException::selfThrow("LTIRelationXML - constructor error : tag " + LTIR_D + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_E)) != NULL)
    ENode = node;
  else
    XMLException::selfThrow("LTIRelationXML - constructor error : tag " + LTIR_E + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_A)) != NULL)
    aNode = node;
  else
    XMLException::selfThrow("LTIRelationXML - constructor error : tag " + LTIR_A + " not found.");
}

LinearTIRXML::~LinearTIRXML()
{
}

void LinearTIRXML::setC(const SiconosMatrix& matrix)
{
  if (CNode == NULL)
  {
    CNode = SiconosDOMTreeTools::createMatrixNode(rootRelationXMLNode, LTIR_C, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(CNode, matrix);
}

void LinearTIRXML::setD(const SiconosMatrix& matrix)
{
  if (DNode == NULL)
  {
    DNode = SiconosDOMTreeTools::createMatrixNode(rootRelationXMLNode, LTIR_D, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(DNode, matrix);
}

void LinearTIRXML::setE(const SiconosMatrix &matrix)
{
  if (ENode == NULL)
  {
    ENode = SiconosDOMTreeTools::createMatrixNode(rootRelationXMLNode, LTIR_E, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(ENode, matrix);
}

void LinearTIRXML::setA(const SiconosVector& vec)
{
  if (aNode == NULL)
  {
    aNode = SiconosDOMTreeTools::createVectorNode(rootRelationXMLNode, LTIR_A, vec);
  }
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(aNode, vec);
}

