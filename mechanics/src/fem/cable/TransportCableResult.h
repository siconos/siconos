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

  int puller12idx{-1};
  int puller21idx{-1};
  Ropeway ropes_up;
  Ropeway ropes_down;

  std::vector<std::shared_ptr<Support>> supports;

  std::vector<Point> q = {};       // positions
  std::vector<Point> R = {};       // internal forces [x,y,z]-> [H,V,B]
  std::vector<double> TS = {};     // tension
  std::vector<int> contacts = {};  // contact points

  int nb_nodes{0};
  double length{0.};
  double elem_length{0.};
  // à convertir en siconos (vecteur ou matrice)
  std::vector<double> punct = {};

  std::vector<double> g = {};

  std::vector<std::vector<Point>> G = {};
  std::vector<std::vector<Point>> T = {};

  /**  */
  std::shared_ptr<siconos::algebra::SiconosVector> q0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> v0{nullptr};

  std::shared_ptr<siconos::algebra::SiconosMatrix> mass{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> fext{nullptr};
};
}  // namespace siconos::fem::cable
