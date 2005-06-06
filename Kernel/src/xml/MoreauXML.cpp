#include "MoreauXML.h"
using namespace std;

MoreauXML::MoreauXML(): OneStepIntegratorXML(), WNode(NULL), ThetaNode(NULL)
{}

MoreauXML::MoreauXML(xmlNode * MoreauNode, map<int, bool> definedDSNumbers):
  OneStepIntegratorXML(MoreauNode, definedDSNumbers), WNode(NULL), ThetaNode(NULL)
{
  xmlNode *node;
  if ((node = SiconosDOMTreeTools::findNodeChild(MoreauNode, MOREAU_W)) != NULL)
    WNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(MoreauNode, MOREAU_THETA)) != NULL)
    ThetaNode = node;
}

MoreauXML::~MoreauXML()
{}

