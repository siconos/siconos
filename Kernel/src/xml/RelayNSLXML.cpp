
#include "RelayNSLXML.h"


RelayNSLXML::RelayNSLXML()
{
  this->CNode = NULL;
  this->DNode = NULL;
}

RelayNSLXML::RelayNSLXML(xmlNode * relayNSLawNode)
  : NonSmoothLawXML(relayNSLawNode)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(relayNSLawNode, RNSL_C)) != NULL)
  {
    this->CNode = node;
  }
  else
  {
    XMLException::selfThrow("RelayNSLawXML - constructor error : tag " + RNSL_C + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(relayNSLawNode, RNSL_D)) != NULL)
  {
    this->DNode = node;
  }
  else
  {
    XMLException::selfThrow("RelayNSLawXML - constructor error : tag " + RNSL_D + " not found.");
  }
}
