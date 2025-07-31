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

/*! \file TransportCableResult.h
  The FEM description of the whole setup (cables, positions ...)
*/
#pragma once

#include <nlohmann/json.hpp>
#include <vector>

#include "Ropeway.h"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::fem::cable {

class Point;

/**  */
class TransportCableResult {
 private:
  TransportCableResult(const TransportCableResult &) = delete;
  TransportCableResult(TransportCableResult &&) = delete;
  TransportCableResult &operator=(const TransportCableResult &) = delete;
  TransportCableResult &operator=(TransportCableResult &&) = delete;

 public:
  /** Default and only constructor */
  TransportCableResult();

  virtual ~TransportCableResult() noexcept = default;

  /** Build supports  */
  void prepareSupport();
  void prepareIneqConstraint(int nb_nodes);

  int exportTC(const std::string &a_fileName, nlohmann::ordered_json &a_output,
               const std::string &a_option = "all");

  int to_json(nlohmann::ordered_json &j, const std::string &a_option = "all");

  /** Vector of all supports: contacts with pylons and stations (pulleys)
  in both up and down ropeways */
  std::vector<std::shared_ptr<Support>> supports;

  /** id of the top pulley in the supports vector */
  int topPulleyId{-1};

  /** id of the down pulley in the supports vector */
  int downPulleyId{-1};

  /** Set of ropes corresponding to the up-way cable */
  Ropeway ropes_up;

  /** Set of ropes corresponding to the down-way cable */
  Ropeway ropes_down;

  /** nodes positions */
  siconos::algebra::SiconosVector q = {};

  /** */
  siconos::algebra::SiconosVector R = {};   // internal forces [x,y,z]-> [H,V,B]
  siconos::algebra::SiconosVector TS = {};  // tension

  /**  */
  std::vector<int> contacts = {};  // contact points

  /** Total number of elements in the finite element mesh */
  int numberOfElements{0};

  /** Total lenght of the cable (wrappers around pulleys + top and down cables) */
  double totalLength{0.};

  /** Length of an element */
  double elementLength{0.};

  /** scalar weight computed at each node (due to carriers) */
  siconos::algebra::SiconosVector weightVector = {};

  /** g(q), vector of the constraints applied to the cable */
  siconos::algebra::SiconosVector gVector = {};

  /** \f$ \nabla_q g(q) \f$ */
  // std::shared_ptr<siconos::algebra::SiconosDenseMatrix> jacobian_g_Over_q{nullptr};
  siconos::algebra::SiconosDenseMatrix jacobian_g_Over_q = {};
  // siconos::algebra::SiconosSparseMatrix jacobian_g_Over_q = {};

  std::shared_ptr<siconos::algebra::SiconosDenseMatrix> T{nullptr};

  /**  */
  std::shared_ptr<siconos::algebra::SiconosVector> q0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> v0{nullptr};

  std::shared_ptr<siconos::algebra::SiconosSparseMatrix> mass{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> fext{nullptr};
};
}  // namespace siconos::fem::cable

// // // Serialization
// // namespace nlohmann {
// // template <>
// // struct adl_serializer<siconos::fem::cable::TransportCableResult> {
// //   static void to_json(json &j, const siconos::fem::cable::TransportCableResult &input) {
// //     j = json{{"EA", input.crossSectionRigidity()},
// //              {"linearDensity", input.linearDensity()},
// //              {"T0", input.initialTension()}};
// //   }
// // };
// // };  // namespace nlohmann
