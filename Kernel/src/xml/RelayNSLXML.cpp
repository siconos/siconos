
#include "RelayNSLXML.h"
using namespace std;

RelayNSLXML::RelayNSLXML(): CNode(NULL), DNode(NULL)
{}

RelayNSLXML::RelayNSLXML(xmlNode * relayNSLawNode)
  : NonSmoothLawXML(relayNSLawNode), CNode(NULL), DNode(NULL)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(relayNSLawNode, RNSL_C)) != NULL)
    CNode = node;
  else
    XMLException::selfThrow("RelayNSLawXML - constructor error : tag " + RNSL_C + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(relayNSLawNode, RNSL_D)) != NULL)
    DNode = node;
  else
    XMLException::selfThrow("RelayNSLawXML - constructor error : tag " + RNSL_D + " not found.");
}
