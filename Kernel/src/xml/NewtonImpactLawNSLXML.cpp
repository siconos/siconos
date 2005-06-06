#include "NewtonImpactLawNSLXML.h"
using namespace std;

NewtonImpactLawNSLXML::NewtonImpactLawNSLXML()
{
  this->ENode = NULL;
}

NewtonImpactLawNSLXML::NewtonImpactLawNSLXML(xmlNode * NewtonImpactLawNSLNode)
  : NonSmoothLawXML(NewtonImpactLawNSLNode)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(NewtonImpactLawNSLNode, NEWTON_E)) != NULL)
  {
    this->ENode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactLawNSLXML - constructor error : tag " + NEWTON_E + " not found.");
  }

}
