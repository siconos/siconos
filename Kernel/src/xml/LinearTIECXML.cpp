#include "LinearTIECXML.h"

LinearTIECXML::LinearTIECXML(): LinearECXML()
{}

LinearTIECXML::LinearTIECXML(xmlNode *node, vector<int> definedDSNumbers)
  : LinearECXML(node, definedDSNumbers)
{}

LinearTIECXML::~LinearTIECXML()
{}

