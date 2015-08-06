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
  /** Basic constructor without plugin builder
     \param [in] q0 SP::SiconosVector initial position of the center of mass.
     \param [in] v0 SP::SiconosVector initial velocity.
     \param [in] mass double& ,the mass.
     \param [in] I SP::SimpleMatrix matrix in R^{3,3}
     \param [in] centerOfMass SP::SiconosVector coordinate of the mass center in the just loaded model
     \param [in] BodyName const std::string& , a string for the body name.
     \param [in] CADFile const std::string& , the cad file.
     \param [in] pluginLib const std::string& , the path to the plugin library.
     \param [in] pluginFct const std::string& , the name of the pluged fonction
   */
  MBTB_Body(SP::SiconosVector q0, SP::SiconosVector v0,
            double& mass,SP::SimpleMatrix I,SP::SiconosVector centerOfMass,
            const std::string& BodyName,  const std::string& CADFile,
            const std::string& pluginLib,  const std::string& pluginFct);


  /** Constructor without plugin builder
     \param [in] q0 SP::SiconosVector initial position of the center of mass.
     \param [in] v0 SP::SiconosVector initial velocity.
     \param [in] mass double& ,the mass.
     \param [in] I SP::SimpleMatrix matrix in R^{3,3}
     \param [in] centerOfMass SP::SiconosVector coordinate of the mass center in the just loaded model
     \param [in] BodyName const std::string& , a string for the body name.
     \param [in] CADFile const std::string& , the cad file.
    */
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
