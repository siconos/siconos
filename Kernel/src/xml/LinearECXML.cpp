#include "LinearECXML.h"
using namespace std;

LinearECXML::LinearECXML(): EqualityConstraintXML()
{}

LinearECXML::LinearECXML(xmlNode *node, vector<int> definedDSNumbers)
  : EqualityConstraintXML(node, definedDSNumbers)
{}

LinearECXML::~LinearECXML()
{}

