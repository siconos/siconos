#include "LagrangianLinearECXML.h"

LagrangianLinearECXML::LagrangianLinearECXML(): EqualityConstraintXML()
{
  this->HNode = NULL;
  this->bNode = NULL;
}

LagrangianLinearECXML::LagrangianLinearECXML(xmlNode *ecNode, vector<int> definedDSNumbers)
  : EqualityConstraintXML(ecNode, definedDSNumbers)
{
  /*    xmlNode *node;

      if ((node=SiconosDOMTreeTools::findNodeChild(LLRelationNode, LLR_H)) !=NULL)
      {
      this->HNode=node;
      }
      else
      {
      XMLException::selfThrow("LLRelationXML - constructor error : tag " + LLR_H + " not found.");
      }

      if ((node=SiconosDOMTreeTools::findNodeChild(LLRelationNode, LLR_B)) !=NULL)
      {
      this->bNode=node;
      }
      else
      {
      XMLException::selfThrow("LLRelationXML - constructor error : tag " + LLR_B + " not found.");
      }
  */
}

LagrangianLinearECXML::~LagrangianLinearECXML()
{}

void LagrangianLinearECXML::setH(SiconosMatrix *matrix)
{
  /*  if( this->HNode == NULL )
    {
      this->HNode = SiconosDOMTreeTools::createMatrixNode(this->rootNode, LLR_H, matrix);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->HNode, matrix);
  */
}

void LagrangianLinearECXML::setB(SiconosVector *vector)
{
  /*  if( this->bNode == NULL )
    {
      this->bNode = SiconosDOMTreeTools::createVectorNode(this->rootNode, LLR_B, vector);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->bNode, vector);
  */
}
