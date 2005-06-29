
#include "TimeDiscretisationXML.h"
using namespace std;


TimeDiscretisationXML::TimeDiscretisationXML():
  rootNode(NULL), hNode(NULL), NNode(NULL), tkNode(NULL),
  hMinNode(NULL), hMaxNode(NULL)
{}

TimeDiscretisationXML::TimeDiscretisationXML(xmlNode * timeDiscretisationNode):
  rootNode(timeDiscretisationNode), hNode(NULL), NNode(NULL), tkNode(NULL),
  hMinNode(NULL), hMaxNode(NULL)
{
  xmlNode *node;
  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_H)) != NULL)
    hNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_N)) != NULL)
    NNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_TK)) != NULL)
    tkNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_HMIN)) != NULL)
    hMinNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_HMAX)) != NULL)
    hMaxNode = node;
}

void TimeDiscretisationXML::updateTimeDiscretisationXML(xmlNode* node, TimeDiscretisation* td)
{
  IN("TimeDiscretisation::updateTimeDiscretisationXML\n");
  rootNode = node;
  OUT("TimeDiscretisation::updateTimeDiscretisationXML\n");
}

