//$Id: NewtonImpactFrictionNSLXML.cpp,v 1.1 2005/03/22 15:55:05 jbarbier Exp $
#include "NewtonImpactFrictionNSLXML.h"


NewtonImpactFrictionNSLXML::NewtonImpactFrictionNSLXML()
{
  this->enNode = NULL;
  this->etNode = NULL;
  this->muNode = NULL;
}

NewtonImpactFrictionNSLXML::NewtonImpactFrictionNSLXML(xmlNode * nslNode)
  : NonSmoothLawXML(nslNode)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, NEWTON_EN)) != NULL)
  {
    this->enNode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag " + NEWTON_EN + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, NEWTON_ET)) != NULL)
  {
    this->etNode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag " + NEWTON_ET + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, NEWTON_MU)) != NULL)
  {
    this->muNode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag " + NEWTON_MU + " not found.");
  }

}
//$Log: NewtonImpactFrictionNSLXML.cpp,v $
//Revision 1.1  2005/03/22 15:55:05  jbarbier
//- class NewtonImpactFriction non smooth law added to the kernel
//
//- xml schema modified for this new class
//- xml schema modified to accept a "joker" for further use of a LMGC90 mechanical plugin
//
//- new test added for the loading/saving of a NewtonImpactFrictionNSL
//
