#include "LinearDSXML.h"
using namespace std;

LinearDSXML::LinearDSXML() :
  DSXML(), ANode(NULL), bNode(NULL), uNode(NULL), ENode(NULL)
{}

LinearDSXML::LinearDSXML(xmlNode * LinearDSNode, const bool& isBVP):
  DSXML(LinearDSNode, isBVP), ANode(NULL), bNode(NULL), uNode(NULL), ENode(NULL)
{
  xmlNode *node;
  // The only required node is A
  node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LDS_A);
  if (node != NULL)
    ANode = node;
  else
    XMLException::selfThrow("LinearDSXML - loadLinearDSProperties error : tag " + LDS_A + " not found.");
  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LDS_B)) != NULL)
    bNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, "uSize")) != NULL)
    uSizeNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LDS_E)) != NULL)
    ENode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LDS_U)) != NULL)
    uNode = node;

}

LinearDSXML::~LinearDSXML()
{}

void LinearDSXML::updateDynamicalSystemXML(xmlNode* rootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  IN("LinearSystemDynamicalSystem::updateDynamicalSystemXML\n");
  ANode = NULL;
  ENode = NULL;
  bNode = NULL;
  uNode = NULL;
  rootDSXMLNode = rootDSXMLNode;
  loadDS(ds);
  OUT("LinearSystemDynamicalSystem::updateDynamicalSystemXML\n");
}

