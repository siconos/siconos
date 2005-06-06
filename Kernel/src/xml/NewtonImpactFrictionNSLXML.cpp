#include "NewtonImpactFrictionNSLXML.h"
using namespace std;

NewtonImpactFrictionNSLXML::NewtonImpactFrictionNSLXML()
{
  this->enNode = NULL;
  this->etNode = NULL;
  this->muNode = NULL;
}

NewtonImpactFrictionNSLXML::NewtonImpactFrictionNSLXML(xmlNode * nslNode)
  : NonSmoothLawXML(nslNode)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, NEWTON_EN)) != NULL)
  {
    this->enNode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag " + NEWTON_EN + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, NEWTON_ET)) != NULL)
  {
    this->etNode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag " + NEWTON_ET + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, NEWTON_MU)) != NULL)
  {
    this->muNode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag " + NEWTON_MU + " not found.");
  }

}
