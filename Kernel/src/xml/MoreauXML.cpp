//$Id: MoreauXML.cpp,v 1.16 2004/09/23 14:09:24 jbarbier Exp $

#include "MoreauXML.h"

MoreauXML::MoreauXML()
{
  this->WNode = NULL;
  this->ThetaNode = NULL;
}

MoreauXML::MoreauXML(xmlNode * MoreauNode, map<int, bool> definedDSNumbers)
  : OneStepIntegratorXML(MoreauNode, definedDSNumbers)
{
  xmlNode *node;
  if ((node = SiconosDOMTreeTools::findNodeChild(MoreauNode, MOREAU_W)) != NULL)
  {
    this->WNode = node;
  }
  else
  {
    this->WNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(MoreauNode, MOREAU_THETA)) != NULL)
  {
    this->ThetaNode = node;
  }
  else
  {
    this->ThetaNode = NULL;
  }
}

MoreauXML::~MoreauXML()
{}

//$Log: MoreauXML.cpp,v $
//Revision 1.16  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.15  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.14  2004/07/12 13:04:34  jbarbier
//- $id$, $log$, $date$, ... added in the XML management files
//- n id calculated with ndof for Lagrangian dynamical systems
//
