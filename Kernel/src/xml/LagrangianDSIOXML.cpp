
#include "LagrangianDSIOXML.h"


LagrangianDSIOXML::LagrangianDSIOXML(): DSInputOutputXML()
{}

LagrangianDSIOXML::LagrangianDSIOXML(xmlNode * dsioNode/*, vector<int> definedDSNumbers */)
  : DSInputOutputXML(dsioNode/*, definedDSNumbers */)
{}

LagrangianDSIOXML::~LagrangianDSIOXML()
{}

