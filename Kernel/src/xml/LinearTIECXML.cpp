//$Id: LinearTIECXML.cpp,v 1.3 2005/01/26 13:50:40 jbarbier Exp $
#include "LinearTIECXML.h"

LinearTIECXML::LinearTIECXML(): LinearECXML()
{}

LinearTIECXML::LinearTIECXML(xmlNode *node, vector<int> definedDSNumbers)
  : LinearECXML(node, definedDSNumbers)
{}

LinearTIECXML::~LinearTIECXML()
{}

//$Log: LinearTIECXML.cpp,v $
//Revision 1.3  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.2  2005/01/17 14:09:34  jbarbier
//- LagrangianECXML class added
//
//Revision 1.1  2005/01/17 10:56:27  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//