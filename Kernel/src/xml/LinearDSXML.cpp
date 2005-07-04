#include "LinearDSXML.h"
using namespace std;

LinearDSXML::LinearDSXML() :
  DSXML(), ANode(NULL), fNode(NULL), uNode(NULL), BNode(NULL)
{}

LinearDSXML::LinearDSXML(xmlNode * LinearDSNode, const bool& isBVP):
  DSXML(LinearDSNode, isBVP), ANode(NULL), fNode(NULL), uNode(NULL), BNode(NULL)
{
  xmlNode *node;
  // The only required node is A
  node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LDS_A);
  if (node != NULL)
    ANode = node;
  else
    XMLException::selfThrow("LinearDSXML - loadLinearDSProperties error : tag " + LDS_A + " not found.");
  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LDS_F)) != NULL)
    fNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, "uSize")) != NULL)
    uSizeNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LDS_B)) != NULL)
    BNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LDS_U)) != NULL)
    uNode = node;

}

LinearDSXML::~LinearDSXML()
{}

void LinearDSXML::updateDynamicalSystemXML(xmlNode* rootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  IN("LinearSystemDynamicalSystem::updateDynamicalSystemXML\n");
  ANode = NULL;
  BNode = NULL;
  fNode = NULL;
  uNode = NULL;
  rootDSXMLNode = rootDSXMLNode;
  loadDS(ds);
  OUT("LinearSystemDynamicalSystem::updateDynamicalSystemXML\n");
}

