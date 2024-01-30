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
/*! \file TransportCableProfil.h

*/

#pragma once

#include <nlohmann/json.hpp>

namespace siconos::fem::cable {

class TransportCableModel;
class TransportCableResult;
class Point;

class TransportCableProfil {
 public:
  TransportCableProfil(const TransportCableModel &a_model, TransportCableResult &a_results);

  virtual ~TransportCableProfil() noexcept = default;

  void computeInitialProfil(int nb_nodes, double a_tol = 1e-20, int a_nmax = 20);

  void computeFEM(int nb_elem, double a_eps = 0.1, double a_tol = 1e-3);

  int to_json(nlohmann::ordered_json &out);

 private:
  const TransportCableModel &r_model;
  TransportCableResult &r_results;

  void compute_punct_load(int nb_elem, double Lc);

  /**
     \param a_X vector of positions, coordinates of cable 'particles'
     \param a_tol tolerance used to activate contacts
   */
  void compute_ineq_constraint(const std::vector<Point> &a_X, double a_tol = 1e-3);
};
}  // namespace siconos::fem::cable
