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
 * limitations under the License.·
 */
#include "Ropeway.h"

#include "MechanicalProperties.h"
#include "PulleyWrapping.h"
#include "Pylon.h"
#include "Rope.h"
#include "SiconosVector.hpp"

void siconos::fem::cable::Ropeway::computeCatenary(const MechanicalProperties &a_meca,
                                                   const std::vector<Pylon> &a_piles,
                                                   int nb_nodes, double a_tol, int a_nmax) {
  /*

  Build ropes and applies Catenary equations for each rope in the ropeway.

  Computes
  - the positions (a 3*nb_nodes vector [x,y,z,...,x,y,z]),
  - the tension (a nb_nodes vector [T,...T])
  - the internal force vector (3*nb_nodes vector [H,V,B,...,H,V,B])
  - the support reaction of piles (6*nb_piles vector [px,py,pz,H,V,B])

   Moreover, it provides with admissibilities for each subspan with a np.array which i-th row
   contains cable unknowns for each subspan np.array([L, etaY, etaZ])
  */
  size_t n = a_piles.size();
  ropes_.clear();
  double T0 = a_meca.initialTension();
  siconos::algebra::SiconosVector3 R0;
  R0.setZero();
  for (size_t k = 0; k < n; k++) {
    size_t k1 = k;
    if (k < n - 1) {  // Last rope seems to be a 'fake' one between last and last point
      k1++;
    }

    ropes_.emplace_back(a_piles[k], a_piles[k1], a_meca, T0, R0, nb_nodes, a_tol, a_nmax);
    // Set T0, R0 for next ropes
    // next T0 = Tension at the last node of the current rope
    T0 = ropes_.back().getTensionAtLastNode();
    // next R0 = R at the last node of the current rope
    if (k < n - 1) R0 = ropes_.back().getLastR();
  }
}

void siconos::fem::cable::Ropeway::prepareSupport(
    std::vector<std::shared_ptr<Support>> &a_supports, int &a_pulleyIdx) const {
  if (!is_down_) {
    for (auto &r : ropes_) {
      addSupport(r, a_supports, a_pulleyIdx);
    }
  } else {
    for (auto r = ropes_.rbegin(); r != ropes_.rend(); r++) {
      addSupport(*r, a_supports, a_pulleyIdx);
    }
  }
}

void siconos::fem::cable::Ropeway::addSupport(
    const Rope &a_rope, std::vector<std::shared_ptr<Support>> &a_supports,
    int &a_pulleyIdx) const {
  if (!a_rope.start_pylon().isStation()) {  // standard case
    a_supports.push_back(std::make_shared<Support>(a_rope.start_pylon().coords(),
                                                   a_rope.start_pylon().get_radius()));
    a_supports.back()->prepare(a_rope);
  } else if (a_pulleyIdx <
             0) {  // station case (--> PulleyWrapping), only if not already in the set
    assert(std::in_range<int>(a_supports.size()));
    a_pulleyIdx = static_cast<int>(a_supports.size());
    a_supports.push_back(std::make_shared<PulleyWrapping>(a_rope.start_pylon().coords()));
  }
}

int siconos::fem::cable::Ropeway::computeNumberOfElements(double element_length,
                                                          double ropewayLength) {
  int nbelem = 0;
  for (auto &r : ropes_) {
    nbelem += r.computeNumberOfElements(element_length, ropewayLength);
  }
  return nbelem;
}

int siconos::fem::cable::Ropeway::initializeFEM(siconos::algebra::SiconosVector &a_q,
                                                siconos::algebra::SiconosVector &a_R,
                                                siconos::algebra::SiconosVector &a_TS,
                                                int q_offset) const {
  int offset = q_offset;
  if (!is_down_) {
    for (auto &r : ropes_) {
      // Compute the local mesh (for each rope) and get back the offset (number of elements in
      // the rope) for the position in q, R and T vectors
      offset += r.initializeFEM(a_q, a_R, a_TS, offset);
    }
  } else {
    for (auto r = ropes_.rbegin(); r != ropes_.rend(); r++) {
      offset += r->initializeFEM(a_q, a_R, a_TS, offset, true);
    }
  }
  return offset;
}

const siconos::fem::cable::Pylon &siconos::fem::cable::Ropeway::getFirstPylon() {
  return ropes_.front().start_pylon();
}

const siconos::fem::cable::Pylon &siconos::fem::cable::Ropeway::getLastPylon() {
  return ropes_.back().start_pylon();
}

double siconos::fem::cable::Ropeway::initialTension() {
  if (ropes_.size()) {
    return ropes_.front().initialTension();
  } else
    return 0.0;
}

double siconos::fem::cable::Ropeway::getTensionAtLastNode() {
  if (ropes_.size()) {
    return ropes_.back().getTensionAtLastNode();
  } else
    return 0.0;
}

double siconos::fem::cable::Ropeway::length() const {
  double l = 0.0;
  // Summation of the lengths of all ropes, from station to station
  // Warning: this does not include the length of the cable around
  // the pulleys at the stations.
  for (auto &r : ropes_) {
    l += r.length();
  }
  return l;
}

const siconos::fem::cable::MechanicalProperties &
siconos::fem::cable::Ropeway::mechanicalProperties0() const {
  return ropes_.front().mechanicalProperties();
}

// int siconos::fem::cable::Ropeway::to_json(nlohmann::ordered_json &j) {
//   // Initialize json array-like fields
//   j = {{"catenaryUnknowns", nlohmann::ordered_json::array()},
//        {"q", nlohmann::ordered_json::array()},
//        {"TS", nlohmann::ordered_json::array()},
//        {"R", nlohmann::ordered_json::array()},
//        {"SR", nlohmann::ordered_json::array()},
//        {"meca_global", nlohmann::ordered_json::array()}};
//   for (auto &r : ropes_) {
//     r.to_json(j);
//   }
//   return 0;
// }

void siconos::fem::cable::Ropeway::set_Down(bool a_value) { is_down_ = a_value; }
