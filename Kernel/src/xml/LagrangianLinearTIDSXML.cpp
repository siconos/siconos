#include "LagrangianLinearTIDSXML.h"
using namespace std;

LagrangianLinearTIDSXML::LagrangianLinearTIDSXML()
  : LagrangianDSXML(), KNode(NULL), CNode(NULL)
{}

LagrangianLinearTIDSXML::LagrangianLinearTIDSXML(xmlNode * LagrangianLinearTIDSNode, bool isBVP):
  LagrangianDSXML(LagrangianLinearTIDSNode, isBVP), KNode(NULL), CNode(NULL)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LTIDS_K)) != NULL)
    KNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, LTIDS_C)) != NULL)
    CNode = node;
}

void LagrangianLinearTIDSXML::updateDynamicalSystemXML(xmlNode* newRootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  IN("LagrangianLinearTIDSXML::updateDynamicalSystemXML\n");
  rootDSXMLNode = newRootDSXMLNode;
  loadDS(ds);
  OUT("LagrangianLinearTIDSXML::updateDynamicalSystemXML\n");
}

