#include "LagrangianLinearECXML.h"

LagrangianLinearECXML::LagrangianLinearECXML(): EqualityConstraintXML()
{}

LagrangianLinearECXML::LagrangianLinearECXML(xmlNode *node, vector<int> definedDSNumbers)
  : EqualityConstraintXML(node, definedDSNumbers)
{}

LagrangianLinearECXML::~LagrangianLinearECXML()
{}

