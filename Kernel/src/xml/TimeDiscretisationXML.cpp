
#include "TimeDiscretisationXML.h"

#include "check.h"


TimeDiscretisationXML::TimeDiscretisationXML()
{
  this->hMin = false;
  this->hMax = false;
  this->hMinNode = false;
  this->hMaxNode = false;
  this->hNode = NULL;
  this->NNode = NULL;
  this->tkNode = NULL;
}

TimeDiscretisationXML::TimeDiscretisationXML(xmlNode * timeDiscretisationNode)
{
  this->loadTimeDiscretisationProperties(timeDiscretisationNode);
}


void TimeDiscretisationXML::loadTimeDiscretisationProperties(xmlNode * timeDiscretisationNode)
{
  xmlNode *node;

  rootNode = timeDiscretisationNode;
  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_H)) != NULL)
  {
    this->hNode = node;
  }
  else
  {
    //XMLException::selfThrow("TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " + TD_H + " not found.");
    cout << "TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " << TD_H << " not found.Optional attribute" << endl;
    this->hNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_N)) != NULL)
  {
    this->NNode = node;
  }
  else
  {
    //XMLException::selfThrow("TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " + TD_N + " not found.");
    cout << "TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " << TD_N << " not found.Optional attribute" << endl;
    this->NNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_TK)) != NULL)
  {
    this->tkNode = node;
  }
  else
  {
    //XMLException::selfThrow("TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " + TD_TK + " not found.");
    cout << "TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " << TD_TK << " not found.Optional attribute" << endl;
    this->tkNode = NULL;
  }

  //Optional attribute
  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_HMIN)) != NULL)
  {
    this->hMinNode = node;
    this->hMin = true;
  }
  else this->hMin = false;
  //Optional attribute
  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_HMAX)) != NULL)
  {
    this->hMaxNode = node;
    this->hMax = true;
  }
  else this->hMax = false;
}

void TimeDiscretisationXML::updateTimeDiscretisationXML(xmlNode* node, TimeDiscretisation* td)
{
  IN("TimeDiscretisation::updateTimeDiscretisationXML\n");
  this->rootNode = node;
  //  this->loadNSDS( nsds );
  OUT("TimeDiscretisation::updateTimeDiscretisationXML\n");
}

