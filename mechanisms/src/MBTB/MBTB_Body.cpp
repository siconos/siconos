#include "MBTB_Body.hpp"

MBTB_Body::MBTB_Body():_mBodyName(""),_cadFileName("")
{
  ;
}
/*
 *
 *
 *
 *
 **/
MBTB_Body::MBTB_Body(SP::SiconosVector q0, SP::SiconosVector v0, double& mass,
                     SP::SimpleMatrix I,SP::SiconosVector centerOfMass,
                     const std::string& BodyName, const  std::string& CADFile,
                     const std::string& pluginLib,  const std::string& pluginFct):
  NewtonEulerDS(q0,v0,mass,I),_mBodyName(BodyName),_cadFileName(CADFile)
{
  setComputeFExtFunction(pluginLib,pluginFct);
  _centerOfMass=centerOfMass;
}

MBTB_Body::MBTB_Body(SP::SiconosVector q0, SP::SiconosVector v0, double& mass,
                     SP::SimpleMatrix I,SP::SiconosVector centerOfMass,
                     const std::string& BodyName, const  std::string& CADFile):
  NewtonEulerDS(q0,v0,mass,I),_mBodyName(BodyName),_cadFileName(CADFile)
{
  _centerOfMass=centerOfMass;

}




MBTB_Body::~MBTB_Body()
{
  ;
}
