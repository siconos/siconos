#include "NonSmoothLawXML.h"
using namespace std;

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

