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
/*! \file TransportCableManager.h

*/

#pragma once
#include <nlohmann/json.hpp>
#include <string>

#include "TransportCableModel.h"
#include "TransportCableResult.h"

/** API to the cable model (TransportCableModel) and its FEM counterpart (TransportCableResult)
 */
namespace siconos::fem::cable {

class TransportCableManager {
 private:
  // Rule of five
  TransportCableManager(const TransportCableManager &) = delete;
  TransportCableManager(TransportCableManager &&) = delete;
  TransportCableManager &operator=(const TransportCableManager &) = delete;
  TransportCableManager &operator=(TransportCableManager &&) = delete;

 public:
  TransportCableManager() = default;
  TransportCableManager(const std::string &a_filename);

  ~TransportCableManager() noexcept = default;

  int importModel(const nlohmann::json &a_input, const std::string &a_filename = "");
  int computeFEM(const nlohmann::json &a_args, const std::string &a_outfile,
                 nlohmann::ordered_json &output);
  int exportTC(const nlohmann::json &a_args, const std::string &a_outfile,
               nlohmann::ordered_json &output);

  int simulation(const nlohmann::json &a_model, const nlohmann::json &a_args,
                 const std::string &a_filename, const std::string &a_outfile,
                 nlohmann::ordered_json &output);

 private:
  TransportCableModel m_model;
  TransportCableResult m_results;

  void computeDS(double a_tolContact = 1e-3, double a_mus = 0.8, double a_mup = 1.1);

  /** Update mass matrix (attribute of m_results)
      \param elem_length elements length
      \param elem_rho linear density
   */
  void compute_mass(double elem_length, double elem_rho);

  /** Update external forces vector (attribute of m_results)
      \param elem_length elements length
      \param elem_rho linear density
  */
  void compute_external_load(double elem_length, double elem_rho);
};
}  // namespace siconos::fem::cable
