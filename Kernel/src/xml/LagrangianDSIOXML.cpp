//$Id: LagrangianDSIOXML.cpp,v 1.3 2005/03/09 15:30:37 jbarbier Exp $

#include "LagrangianDSIOXML.h"


LagrangianDSIOXML::LagrangianDSIOXML(): DSInputOutputXML()
{}

LagrangianDSIOXML::LagrangianDSIOXML(xmlNode * dsioNode/*, vector<int> definedDSNumbers */)
  : DSInputOutputXML(dsioNode/*, definedDSNumbers */)
{}

LagrangianDSIOXML::~LagrangianDSIOXML()
{}

//$Log: LagrangianDSIOXML.cpp,v $
//Revision 1.3  2005/03/09 15:30:37  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.2  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.1  2005/01/17 10:56:27  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//
