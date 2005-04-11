
#include "LinearDSIOXML.h"

LinearDSIOXML::LinearDSIOXML(): DSInputOutputXML()
{
}

LinearDSIOXML::LinearDSIOXML(xmlNode * dsioNode): DSInputOutputXML(dsioNode/*, definedDSNumbers */)
{
  xmlNode *node;

  this->rootDSIOXMLNode = dsioNode;
  if ((node = SiconosDOMTreeTools::findNodeChild(dsioNode, LINEARDSIO_A)) != NULL)
    this->ANode = node;
  else
    XMLException::selfThrow("LinearDSIOXML - constructor error : tag " + LINEARDSIO_A + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(dsioNode, LINEARDSIO_B)) != NULL)
    this->BNode = node;
  else
    XMLException::selfThrow("LinearDSIOXML - constructor error : tag " + LINEARDSIO_B + " not found.");
}

LinearDSIOXML::~LinearDSIOXML()
{}

void LinearDSIOXML::setA(SiconosMatrix *matrix)
{
  if (this->ANode == NULL)
  {
    this->ANode = SiconosDOMTreeTools::createMatrixNode(this->rootDSIOXMLNode, LINEARDSIO_A, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixValue(ANode, matrix);
}

void LinearDSIOXML::setB(SiconosMatrix *matrix)
{
  if (this->BNode == NULL)
  {
    this->BNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSIOXMLNode, LINEARDSIO_B, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixValue(BNode, matrix);
}

//void LinearDSIOXML::setE(SiconosMatrix *matrix)
//{
//  if( this->ENode == NULL )
//  {
//    this->ENode = SiconosDOMTreeTools::createMatrixNode(this->rootDSIOXMLNode, LINEARDSIO_E, matrix);
//  }
//  else SiconosDOMTreeTools::setSiconosMatrixValue(ENode, matrix);
//}
//
//void LinearDSIOXML::setA(SiconosVector *vector)
//{
//  if( this->aNode == NULL )
//  {
//    this->aNode = SiconosDOMTreeTools::createVectorNode(this->rootDSIOXMLNode, LINEARDSIO_A, vector);
//  }
//  else SiconosDOMTreeTools::setSiconosVectorValue(aNode, vector);
//}

