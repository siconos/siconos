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

/*! \file PulleyWrapping.h

  PulleyWrapping class (top and down pylons, at station)
*/

#pragma once

#include "SiconosVector.hpp"
#include "Support.h"

namespace siconos::fem::cable {

/** Description of the cable part around the top and down pulleys, at the stations.
 *  We consider a winding angle of \f$ \pi \f$
 */
class PulleyWrapping : public Support {
 private:
  /** Tension */
  double tension_{0.};

  /** Initial value for theta angle used to compute discretisation of the rope around the
   * pulley */
  double theta_start_{0.};

 public:
  explicit PulleyWrapping(const siconos::algebra::SiconosVector3 &coordinates)
      : Support(coordinates, 0.) {}  // Radius computed later

  virtual ~PulleyWrapping() noexcept = default;

  /** \return the length (after computation) of the cable around the pulley
   *    (considering a winding angle of \f$ \pi \f$)
   */
  double length(const class Ropeway &a_rope) const;

  double tension() const { return tension_; }

  //------------ statique -------------
  virtual void prepare(const Rope &a_rope) override {};
  virtual void prepare(const Pylon &a_start, const Pylon &a_end, double T) override;

  /**
   * @brief Compute nodes coordinates of the discretized rope around the pulley
   *
   * @param[in,out] nb number of nodes used in the discretization
   * @param[in, out] nodes_coords global nodes coordinates vector
   * @param[in] q_offset position (node number, not dof!) in the global vector to be modified
   * @return int
   */
  int computeMesh(int nb, siconos::algebra::SiconosVector &a_q, int q_offset = 0) const;

  virtual void compute(const siconos::algebra::SiconosVector3 &a_p, double a_tol, double &g,
                       Eigen::Ref<siconos::algebra::SiconosVector3> G,
                       Eigen::Ref<siconos::algebra::SiconosVector3> T, int &c) override;
  //------------ dynamique -------------
  virtual bool isContact(const Eigen::Ref<siconos::algebra::SiconosVector3> &a_p,
                         double a_tol) override;

  /** display object to screen */
  void display() const override;
};

}  // namespace siconos::fem::cable

// Serialization
namespace nlohmann {
template <>
struct adl_serializer<siconos::fem::cable::PulleyWrapping> {
  // No need of a from_json for pulleys
  //   static siconos::fem::cable::PulleyWrapping from_json(const json &j) {
  //     return siconos::fem::cable::PulleyWrapping(j);
  //   }

  static void to_json(ordered_json &j, const siconos::fem::cable::PulleyWrapping &input) {
    // Serialize base class part
    adl_serializer<siconos::fem::cable::Support>::to_json(j, input);
    j["tension"] = input.tension();
    // j = json{
    //     {"radius", input.radius()}, {"center", input.center()}, {"tension",
    //     input.tension()}};
  }
};
}  // namespace nlohmann