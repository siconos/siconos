
#include "LinearSystemDSXML.h"

#include "check.h"


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

  node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LSDS_A);
  if (node != NULL)
  {
    this->ANode = node;
  }
  else
  {
    XMLException::selfThrow("LinearSystemDSXML - loadLinearSystemDSProperties error : tag " + LSDS_A + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LSDS_B)) != NULL)
  {
    this->BNode = node;
  }
  else
  {
    XMLException::selfThrow("LinearSystemDSXML - loadLinearSystemDSProperties error : tag " + LSDS_B + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LSDS_U)) != NULL)
  {
    this->uNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LSDS_U + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSXMLNode, LSDS_F)) != NULL)
  {
    this->fNode = node;
  }
  else
  {
    XMLException::selfThrow("LagrangianNLDSXML - loadLagrangianNLDSProperties error : tag " + LSDS_F + " not found.");
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

//$Log: LinearSystemDSXML.cpp,v $
//Revision 1.14  2004/09/10 08:04:50  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.13  2004/08/23 14:30:02  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.12  2004/08/20 15:26:45  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NSDS and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.11  2004/07/29 14:25:43  jbarbier
