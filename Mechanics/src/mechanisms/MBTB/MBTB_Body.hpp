#ifndef MBTB_BODY
#define MBTB_BODY
#include "SiconosKernel.hpp"
/**
 * \brief This class implements a body in a multi-bodies system.
 * It inherits from Siconos::NewtonEulerDS.
 */
class MBTB_Body : public NewtonEulerDS
{
protected:
  /** coordinate of the center of mass in the just loaded model*/
  SP::SiconosVector _centerOfMass;
  //! The name of the body.
  const std::string& _mBodyName;
  //! The cad file.
  const std::string& _cadFileName;
  MBTB_Body();
public:
  /**
     builder
     \param [in] SP::SiconosVector q0, initial position of the center of mass.
     \param [in] SP::SiconosVector v0, initial velocity.
     \param [in]  double& mass,the mass.
     \param [in]  SP::SimpleMatrix I,matrix in R^{3,3}
     \param [in]  SP::SiconosVector centerOfMass,coordinate of the mass center in the just loaded model
     \param [in] const std::string& BodyName, a string for the body name.
     \param [in] const std::string& CADFile, the cad file.
     \param [in] const std::string& pluginLib, the path to the plugin library.
     \param [in] const std::string& plunginFct, the name of the pluged fonction
   */
  MBTB_Body(SP::SiconosVector q0, SP::SiconosVector v0,
            double& mass,SP::SimpleMatrix I,SP::SiconosVector centerOfMass,
            const std::string& BodyName,  const std::string& CADFile,
            const std::string& pluginLib,  const std::string& plunginFct);
  MBTB_Body(SP::SiconosVector q0, SP::SiconosVector v0,
            double& mass,SP::SimpleMatrix I,SP::SiconosVector centerOfMass,
            const std::string& BodyName,  const std::string& CADFile);
  virtual ~MBTB_Body();
  inline SP::SiconosVector centerOfMass() const
  {
    return _centerOfMass;
  }
};
TYPEDEF_SPTR(MBTB_Body);

#endif
