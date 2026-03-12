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

/*! \file Ropeway.h

  Ropeway class
*/

#pragma once

#include <nlohmann/json.hpp>
#include <vector>

#include "Rope.h"
#include "SiconosVector.hpp"

namespace siconos::fem::cable {

class MechanicalProperties;
class Pylon;
class Rope;
class Support;
class Point;

/**  A sequence of ropes used to described the whole cable between two pulleys
 *   (up and down  stations)
 *
 */
class Ropeway {
 private:
  /** All pieces of rope */
  std::vector<Rope> ropes_ = {};

  /** true if the ropeway goes down */
  bool is_down_{false};

  Ropeway(const Ropeway &) = delete;
  Ropeway &operator=(const Ropeway &) = delete;
  Ropeway(Ropeway &&) = delete;
  Ropeway &operator=(Ropeway &&) = delete;

  /** Create new support from the lower pile of a given rope

     \param[in] a_rope the considered Rope
     \param[in,out] a_supports the vector of all supports
     \param[in,out] a_pulleyIdx current number of added support (internal counter)
   */
  void addSupport(const Rope &a_rope, std::vector<std::shared_ptr<Support>> &a_supports,
                  int &a_pulleyIdx) const;

 public:
  /** Default and only constructor */
  Ropeway() = default;

  virtual ~Ropeway() noexcept = default;

  /**  Build the vector of ropes and computes their initial profiles using Catenary equations

     \param a_meca geometric and material properties of the cable
     \param a_piles the vector of pylons, must be ordered!
     \param nb_nodes number of nodes by segment (between two pylons)
     \param a_tol tolerance used to compute admissibility conditions
     \param a_nmax max number of iterations used to compute admissibility conditions
  */
  void computeCatenary(const MechanicalProperties &a_meca, const std::vector<Pylon> &a_piles,
                       int nb_nodes, double a_tol = 1e-20, int a_nmax = 20);

  /**
     Create and prepare all supports, from the list of declared ropes
     (a support for each 'lower' pile of each rope).

     \param[in,out] a_supports the vector of all supports
     \param[in,out] a_pulleyIdx current number of added support (internal counter)
   */
  void prepareSupport(std::vector<std::shared_ptr<Support>> &a_supports,
                      int &a_pulleyIdx) const;

  /** Compute the number of elements to be used in the current cable

      \param[in] element_length expected element size
      \param[in] total_length total length of the cable (all ropes, the whole setup)
  */
  int computeNumberOfElements(double element_length, double total_length);

  /** Computes the profile of each rope in the ropeway
   * \param[out] a_q nodes coordinates
   * \param[out] a_R reaction forces at nodes
   * \param[out] a_TS tension at nodes
   * \param[in] q_offset offset (index) in a_q vector
   */
  int initializeFEM(siconos::algebra::SiconosVector &a_q, siconos::algebra::SiconosVector &a_R,
                    siconos::algebra::SiconosVector &a_TS, int q_offset) const;

  const Pylon &getFirstPylon();
  const Pylon &getLastPylon();
  double initialTension();
  double getTensionAtLastNode();

  /** \return the total length of the ropeway (from station to station but excluding the length
       of the cable around the pulleys at station)
   */
  double length() const;

  const MechanicalProperties &mechanicalProperties0() const;

  void set_Down(bool a_value);

  // For serialisation
  friend struct nlohmann::adl_serializer<Ropeway>;
};

}  // namespace siconos::fem::cable

namespace nlohmann {
template <>
struct adl_serializer<siconos::fem::cable::Ropeway> {
  static siconos::fem::cable::Ropeway from_json(const json &j) {
    return siconos::fem::cable::Ropeway(j);
  }

  static void to_json(ordered_json &j, const siconos::fem::cable::Ropeway &input) {
    j = {{"catenaryUnknowns", nlohmann::ordered_json::array()},
         {"q", nlohmann::ordered_json::array()},
         {"TS", nlohmann::ordered_json::array()},
         {"R", nlohmann::ordered_json::array()},
         {"SR", nlohmann::ordered_json::array()},
         {"meca_global", nlohmann::ordered_json::array()}};
    for (auto &r : input.ropes_) {
      nlohmann::adl_serializer<siconos::fem::cable::Rope>::to_json(j, r);
    }
  }
};
}  // namespace nlohmann