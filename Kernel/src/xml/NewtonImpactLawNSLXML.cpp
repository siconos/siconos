//$Id: NewtonImpactLawNSLXML.cpp,v 1.2 2004/09/21 11:49:10 jbarbier Exp $
#include "NewtonImpactLawNSLXML.h"


NewtonImpactLawNSLXML::NewtonImpactLawNSLXML()
{
  this->ENode = NULL;
}

NewtonImpactLawNSLXML::NewtonImpactLawNSLXML(xmlNode * NewtonImpactLawNSLNode)
  : NonSmoothLawXML(NewtonImpactLawNSLNode)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(NewtonImpactLawNSLNode, NEWTON_E)) != NULL)
  {
    this->ENode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactLawNSLXML - constructor error : tag " + NEWTON_E + " not found.");
  }

}
//$Log: NewtonImpactLawNSLXML.cpp,v $
//Revision 1.2  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.1  2004/07/06 14:54:49  acary
//Renaming NSLaw into NonSmoothLaw
//Renaming RelayNSLaw into RelayNSL
//Renaming CCNSLaw into ComplementarityConditionNSL
//Renaming NewtonImpactLaw into NewtonImpactLawNSL
//
//Revision 1.1  2004/06/30 09:44:35  acary
//Added NewtonImpactLawNSL
//