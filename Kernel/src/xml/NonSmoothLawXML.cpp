//$Id: NonSmoothLawXML.cpp,v 1.5 2004/09/14 13:49:58 jbarbier Exp $

#include "NonSmoothLawXML.h"


NonSmoothLawXML::NonSmoothLawXML()
{}

NonSmoothLawXML::NonSmoothLawXML(xmlNode *node)
{
  this->rootNSLawXMLNode = node;
}

NonSmoothLawXML::~NonSmoothLawXML()
{}

void NonSmoothLawXML::updateNonSmoothLawXML(xmlNode* node, NonSmoothLaw* nsl)
{
  IN("NonSmoothLawXML::updateNonSmoothLawXML\n");
  this->rootNSLawXMLNode = node;
  OUT("RelaNonSmoothLawXMLtionXML::updateNonSmoothLawXML\n");
}

//$Log: NonSmoothLawXML.cpp,v $
//Revision 1.5  2004/09/14 13:49:58  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.4  2004/07/29 14:25:44  jbarbier
//- $Log: NonSmoothLawXML.cpp,v $
//- Revision 1.5  2004/09/14 13:49:58  jbarbier
//- - files added in sample/ to run run the main_siconos test program
//-
//- - all the platform can now be saved in an XML file when it is created manually
//- and $Id: NonSmoothLawXML.cpp,v 1.5 2004/09/14 13:49:58 jbarbier Exp $ added
//
