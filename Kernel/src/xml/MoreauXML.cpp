
#include "MoreauXML.h"

MoreauXML::MoreauXML()
{
  this->WNode = NULL;
  this->ThetaNode = NULL;
}

MoreauXML::MoreauXML(xmlNode * MoreauNode, map<int, bool> definedDSNumbers)
  : OneStepIntegratorXML(MoreauNode, definedDSNumbers)
{
  xmlNode *node;
  if ((node = SiconosDOMTreeTools::findNodeChild(MoreauNode, MOREAU_W)) != NULL)
  {
    this->WNode = node;
  }
  else
  {
    this->WNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(MoreauNode, MOREAU_THETA)) != NULL)
  {
    this->ThetaNode = node;
  }
  else
  {
    this->ThetaNode = NULL;
  }
}

MoreauXML::~MoreauXML()
{}

