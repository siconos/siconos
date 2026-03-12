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

/*! \file Rope.h

  Rope class
*/
#pragma once

#include "MechanicalProperties.h"
#include "SiconosVector.hpp"

namespace siconos::fem::cable {

class Pylon;

/** Describes a  discrete piece of cable, between two pylons.
 *
 */
class Rope {
 private:
  /** Unknowns in the catenary equation:  \f$ [L, \eta_y, \eta_z] \f$
        L: unstretched length of the cable
        etaY : initial y-slope (at first node)
        etaZ : initial z-slope (at first node)
 */
  siconos::algebra::SiconosVector3 catenaryUnknowns_;

  /** Coordinates of rope nodes */
  siconos::algebra::SiconosVector nodes_coords_ = {};

  /** Internal forces (at each node) [x,y,z]-> [H,V,B] */
  siconos::algebra::SiconosVector internal_forces_ = {};

  /** Tension (at each node)*/
  siconos::algebra::SiconosVector TS_ = {};

  /** Mechanical properties of the rope (EA, rho ...) */
  MechanicalProperties mechanicalProp_;

  /** First pylon supporting the rope */
  const Pylon &start_pylon_;

  /** Second pylon supporting the rope */
  const Pylon &end_pylon_;

  /** Support reaction at contact [H,V,B]*/
  siconos::algebra::SiconosVector3 support_reaction_ = {0., 0., 0.};

  /** True if the rope is the last one in the line */
  bool m_last_{false};

  /** level of tolerance used to compute admissibility conditions */
  double tol_ = 1e-7;

  /** max number of iterations used in compute admissibility conditions */
  int max_iter_ = 50;

  /** number of elements used to discretize the rope */
  int number_of_elements_{0};

  Rope() = delete;
  Rope &operator=(const Rope &) = delete;
  Rope &operator=(Rope &&) = delete;
  Rope(const Rope &) = delete;

 public:
  /** Build a rope span (a piece of cable between two pylons)
      (compute initial profile of the Rope using Catenary equations)

      Update positions, internal forces and tension in the rope.

      \param start_pylon first pylon supporting the rope
      \param end_pylon second pylon
      \param meca_prop mechanical properties of the rope span
      \param T0 initial tension applied at the beginning of the rope
      \param R0 support reaction at the beginning of the rope
      \param nb_nodes number of nodes uses to compute profile with catenary equations
      \param tol tolerance used to perform (stop) Newton-Raphson iterations
      \param max max. number of iterations used in Newton-Raphson
  */
  explicit Rope(const Pylon &start_pylon, const Pylon &end_pylon,
                MechanicalProperties meca_prop, double T0,
                const siconos::algebra::SiconosVector3 &R0, int nb_nodes, double tol,
                int max_iter);

  Rope(Rope &&) noexcept = default;  // Required for std::vector<Rope> in Ropeway

  ~Rope() noexcept = default;

  /** Compute the number of elements to be used in the rope

      \param[in] element_length expected element size
      \param[in] total_length total length of the cable (the whole setup)
  */
  int computeNumberOfElements(double element_length, double total_length);

  /** Compute a mesh for the current piece of rope
   * \param[out] a_q nodes coordinates
   * \param[out] a_R reaction forces at nodes
   * \param[out] a_TS tension at nodes
   * \param[in] q_offset offset (index) in a_q vector
   * \param[in] a_reverse true when 'down' ropeway is concerned
   * \return the number of elements to be used in the rope
   */
  int initializeFEM(siconos::algebra::SiconosVector &a_q, siconos::algebra::SiconosVector &a_R,
                    siconos::algebra::SiconosVector &a_TS, int q_offset,
                    bool a_reverse = false) const;

  /** \return tension at the first node of the Rope */
  double initialTension();

  /** \return tension at the last node of the Rope */
  double getTensionAtLastNode();

  /** \return reaction at the last node of the Rope */
  Eigen::Ref<const siconos::algebra::SiconosVector3> getLastR() const;

  /** \return the length of the rope span */
  double length() const { return catenaryUnknowns_(0); };

  Eigen::Ref<const siconos::algebra::SiconosVector3> supportReaction() const {
    return support_reaction_;
  };
  const MechanicalProperties &mechanicalProperties() const { return mechanicalProp_; }
  const Pylon &start_pylon() const noexcept { return start_pylon_; };

  /** Print Rope params to screen */
  void display() const;

  // For serialisation
  friend struct nlohmann::adl_serializer<Rope>;
};

// Free functions
/**
   \returns an initial guess for the initial length and the slopes of the rope, between two
   given pylons
   \param[in] start coordinates of the first node in the rope
   \param[in] end coordinates of the last node in the rope
 */
siconos::algebra::SiconosVector3 guess(const siconos::algebra::SiconosVector3 &start,
                                       const siconos::algebra::SiconosVector3 &end);

/** Computes initial lengths and slopes by solving admissibility conditions for the catenary
    See Ch. Bertrand Phd

  \param[in] a_meca geometric and material description of the cable
  \param[in] start coordinates of the first node in the rope
  \param[in] end coordinates of the last node in the rope
  \param max_iter max number of iterations
  \param tol required tolerance
  \return [Unstretched length, initial y-slope, initial z-slope]

*/
siconos::algebra::SiconosVector3 computeAdmissibilityConditions(
    const MechanicalProperties &a_meca, const siconos::algebra::SiconosVector3 &start,
    const siconos::algebra::SiconosVector3 &end, int max_iter, double tol);

/** Discretize the rope by applying catenary equation

  \param[in] a_meca geometric and material description of the cable
  \param[in] end coordinates of the last point in the rope (end pylon)
  \param[in] cable_inc [L, etaY, etaZ], L is the unstretched length of the cable,
    etaY the initial y-slope and etaZ the initial z-slope
  \param[in] nb_nodes number of nodes used to compute catenary
  \param[in, out] computed nodes coordinates, size = (nb_nodes - 1) * 3 [N0x, N0y, N0z,
      N1x, ...], we do not keep last point
  \param[in, out] computed internal forces at each node, size = nb_nodes * 3
      [H, V, B, ..., H, V, B ...]
  \param[in, out] computed tension at each node, size = nb_nodes
  \param[in] q_offset value of the starting index (node number, not dof!) where results are
     computed in positions, internal_forces and tension
  \param[in] a_reverse true to compute catenary starting from the last pylon

*/
void computeCatenary(const MechanicalProperties &a_meca,
                     const siconos::algebra::SiconosVector3 &end,
                     const siconos::algebra::SiconosVector3 &cable_inc, int nb_nodes,
                     siconos::algebra::SiconosVector &positions,
                     siconos::algebra::SiconosVector &internal_forces,
                     siconos::algebra::SiconosVector &tension, int q_offset = 0,
                     bool a_reverse = false);

}  // namespace siconos::fem::cable

// Serialization
namespace nlohmann {
template <>
struct adl_serializer<siconos::fem::cable::Rope> {
  static siconos::fem::cable::Rope from_json(const json &j) {
    return siconos::fem::cable::Rope(j);
  }

  static void to_json(ordered_json &j, const siconos::fem::cable::Rope &input) {
    for (int i = 0; i < input.nodes_coords_.size(); ++i) {
      j["q"].push_back(input.nodes_coords_(i));
    }
    j["SR"].push_back(input.supportReaction());

    if (input.TS_.size()) {
      j["catenaryUnknowns"].push_back(input.catenaryUnknowns_);

      for (int i = 0; i < input.TS_.size(); ++i) {
        j["TS"].push_back(input.TS_(i));
      }
      for (int i = 0; i < input.internal_forces_.size(); ++i) {
        j["R"].push_back(input.internal_forces_(i));
      }

      std::vector<double> buff = {input.mechanicalProperties().initialTension(),
                                  input.mechanicalProperties().crossSectionRigidity(),
                                  input.mechanicalProperties().linearDensity()};
      j["meca_global"].push_back(buff);
    }
  }
};
}  // namespace nlohmann