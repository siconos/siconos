
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
//$Log: RelayNSLXML.cpp,v $
//Revision 1.5  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.4  2004/07/29 14:25:45  jbarbier
