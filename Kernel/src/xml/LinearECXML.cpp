#include "LinearECXML.h"

LinearECXML::LinearECXML(): EqualityConstraintXML()
{}

LinearECXML::LinearECXML(xmlNode *node, vector<int> definedDSNumbers)
  : EqualityConstraintXML(node, definedDSNumbers)
{}

LinearECXML::~LinearECXML()
{}

