#include "LagrangianECXML.h"

LagrangianECXML::LagrangianECXML(): EqualityConstraintXML()
{}

LagrangianECXML::LagrangianECXML(xmlNode *node, vector<int> definedDSNumbers)
  : EqualityConstraintXML(node, definedDSNumbers)
{}

LagrangianECXML::~LagrangianECXML()
{}

