#include "LagrangianECXML.h"
using namespace std;

LagrangianECXML::LagrangianECXML(): EqualityConstraintXML()
{}

LagrangianECXML::LagrangianECXML(xmlNode *node, vector<int> definedDSNumbers)
  : EqualityConstraintXML(node, definedDSNumbers)
{}

LagrangianECXML::~LagrangianECXML()
{}

