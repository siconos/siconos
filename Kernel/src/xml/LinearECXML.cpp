//$Id: LinearECXML.cpp,v 1.2 2005/01/26 13:50:40 jbarbier Exp $
#include "LinearECXML.h"

LinearECXML::LinearECXML(): EqualityConstraintXML()
{}

LinearECXML::LinearECXML(xmlNode *node, vector<int> definedDSNumbers)
  : EqualityConstraintXML(node, definedDSNumbers)
{}

LinearECXML::~LinearECXML()
{}

//$Log: LinearECXML.cpp,v $
//Revision 1.2  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.1  2005/01/17 10:56:27  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//