#include "LinearSystemDSXML.h"
using namespace std;

LinearSystemDSXML::LinearSystemDSXML() : DSXML()
{
  this->ANode = NULL;
  this->BNode = NULL;
  this->fNode = NULL;
  this->uNode = NULL;
}

LinearSystemDSXML::LinearSystemDSXML(xmlNode * LinearSystemDSNode, bool isBVP)
  : DSXML(LinearSystemDSNode, isBVP)
{
  this->ANode = NULL;
  this->BNode = NULL;
  this->fNode = NULL;
  this->uNode = NULL;
  this->loadLinearSystemDSProperties();
}

LinearSystemDSXML::~LinearSystemDSXML()
{}

void LinearSystemDSXML::loadLinearSystemDSProperties()
{
  xmlNode *node;

  node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LDS_A);
  if (node != NULL)
  {
    this->ANode = node;
  }
  else
  {
    XMLException::selfThrow("LinearSystemDSXML - loadLinearSystemDSProperties error : tag " + LDS_A + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LDS_B)) != NULL)
  {
    this->BNode = node;
  }
  else
  {
    XMLException::selfThrow("LinearSystemDSXML - loadLinearSystemDSProperties error : tag " + LDS_B + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LDS_U)) != NULL)
  {
    this->uNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LDS_U + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LDS_F)) != NULL)
  {
    this->fNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianDSXML - loadLagrangianDSProperties error : tag " + LDS_F + " not found.");
  }
}

void LinearSystemDSXML::updateDynamicalSystemXML(xmlNode* rootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  IN("LinearSystemDynamicalSystem::updateDynamicalSystemXML\n");
  this->ANode = NULL;
  this->BNode = NULL;
  this->fNode = NULL;
  this->uNode = NULL;
  this->rootDSXMLNode = rootDSXMLNode;
  this->loadDS(ds);
  OUT("LinearSystemDynamicalSystem::updateDynamicalSystemXML\n");
}

