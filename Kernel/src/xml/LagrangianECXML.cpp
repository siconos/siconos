//$Id: LagrangianECXML.cpp,v 1.2 2005/01/26 13:50:40 jbarbier Exp $
#include "LagrangianECXML.h"

LagrangianECXML::LagrangianECXML(): EqualityConstraintXML()
{}

LagrangianECXML::LagrangianECXML(xmlNode *node, vector<int> definedDSNumbers)
  : EqualityConstraintXML(node, definedDSNumbers)
{}

LagrangianECXML::~LagrangianECXML()
{}

//$Log: LagrangianECXML.cpp,v $
//Revision 1.2  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.1  2005/01/17 14:09:34  jbarbier
//- LagrangianECXML class added
//