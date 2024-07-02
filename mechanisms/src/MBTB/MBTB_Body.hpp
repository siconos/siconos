/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef MBTB_BODY
#define MBTB_BODY

#include "NewtonEulerDS.hpp"  // base class

namespace siconos::mechanisms {

/**
 * \brief This class implements a body in a multi-bodies system.
 * It inherits from Siconos::NewtonEulerDS.
 */
class MBTB_Body : public siconos::modeling::NewtonEulerDS {
 protected:
  /** coordinate of the center of mass in the just loaded model*/
  std::shared_ptr<siconos::algebra::SiconosVector> _centerOfMass{nullptr};
  //! The name of the body.
  const std::string _mBodyName{""};
  //! The cad file.
  const std::string _cadFileName{""};

  MBTB_Body() = default;

 public:
  /** Constructor without plugin builder

      \param [in] q0  initial position of the center of mass.
      \param [in] v0  initial velocity.
      \param [in] mass double& ,the mass.
      \param [in] I  matrix in R^{3,3}
      \param [in] centerOfMass  coordinate of the mass center in
      the just loaded model \param [in] BodyName const std::string& , a string
      for the body name. \param [in] CADFile const std::string& , the cad file.
    */
  MBTB_Body(std::shared_ptr<siconos::algebra::SiconosVector> q0,
            std::shared_ptr<siconos::algebra::SiconosVector> v0, double mass,
            std::shared_ptr<siconos::algebra::SiconosMatrix> I,
            std::shared_ptr<siconos::algebra::SiconosVector> centerOfMass,
            const std::string& BodyName, const std::string& CADFile);

  /** Constructor with plugins

      \param [in] q0  initial position of the center of mass.
      \param [in] v0  initial velocity.
      \param [in] mass double& ,the mass.
      \param [in] I  matrix in R^{3,3}
      \param [in] centerOfMass  coordinate of the mass center in
      the just loaded model \param [in] BodyName const std::string& , a string
      for the body name. \param [in] CADFile const std::string& , the cad file.
      \param [in] pluginLib const std::string& , the path to the plugin library.
      \param [in] pluginFct const std::string& , the name of the pluged fonction
   */
  MBTB_Body(std::shared_ptr<siconos::algebra::SiconosVector> q0,
            std::shared_ptr<siconos::algebra::SiconosVector> v0, double mass,
            std::shared_ptr<siconos::algebra::SiconosMatrix> I,
            std::shared_ptr<siconos::algebra::SiconosVector> centerOfMass,
            const std::string& BodyName, const std::string& CADFile,
            const std::string& pluginLib, const std::string& pluginFct);

  virtual ~MBTB_Body() noexcept = default;

  inline auto centerOfMass() const { return _centerOfMass; }
};

}  // namespace siconos::mechanisms

#endif
