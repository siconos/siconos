#include "LagrangianLinearECXML.h"
using namespace std;

LagrangianLinearECXML::LagrangianLinearECXML(): EqualityConstraintXML()
{
  this->HNode = NULL;
  this->bNode = NULL;
}

LagrangianLinearECXML::LagrangianLinearECXML(xmlNode *ecNode, vector<int> definedDSNumbers)
  : EqualityConstraintXML(ecNode, definedDSNumbers)
{
  /*    xmlNode *node;

      if ((node=SiconosDOMTreeTools::findNodeChild(LLRelationNode, LLEC_H)) !=NULL)
      {
      this->HNode=node;
      }
      else
      {
      XMLException::selfThrow("LLRelationXML - constructor error : tag " + LLEC_H + " not found.");
      }

      if ((node=SiconosDOMTreeTools::findNodeChild(LLRelationNode, LLEC_B)) !=NULL)
      {
      this->bNode=node;
      }
      else
      {
      XMLException::selfThrow("LLRelationXML - constructor error : tag " + LLEC_B + " not found.");
      }
  */
}

LagrangianLinearECXML::~LagrangianLinearECXML()
{}

void LagrangianLinearECXML::setH(SiconosMatrix *matrix)
{
  /*  if( this->HNode == NULL )
    {
      this->HNode = SiconosDOMTreeTools::createMatrixNode(this->rootNode, LLEC_H, matrix);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->HNode, matrix);
  */
}

void LagrangianLinearECXML::setB(SiconosVector *vector)
{
  /*  if( this->bNode == NULL )
    {
      this->bNode = SiconosDOMTreeTools::createVectorNode(this->rootNode, LLEC_B, vector);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->bNode, vector);
  */
}
