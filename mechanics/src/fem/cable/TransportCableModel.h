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

/*! \file TransportCableModel.h
  Full description of the setup (cable, piles, carriers)
  before discretization
*/
#pragma once

#include <vector>

#include "Cable.h"
#include "Carriers.h"
#include "Pile.h"

namespace siconos::fem::cable {

class TransportCableModel {
 public:
  TransportCableModel() = default;

  /** Read data from a json file
      \param a_filename input file name (json)
  */
  TransportCableModel(const std::string &a_filename);

  virtual ~TransportCableModel() noexcept = default;

  int from_file(const std::string &a_fileName);
  int from_json(const nlohmann::json &j);
  bool isLoaded();

  const Cable &get_cable() const;
  const Carriers &get_carriers() const;
  const std::vector<Pile> &get_piles1() const;
  const std::vector<Pile> &get_piles2() const;

 private:
  Cable m_cable;
  Carriers m_carriers;
  Pile m_stationUp;  // drive station
  Pile m_stationDown;
  std::vector<Pile> m_piles = {};  // Roller batteries
  // Rule of five
  TransportCableModel(const TransportCableModel &) = delete;
  TransportCableModel(TransportCableModel &&) = delete;
  TransportCableModel &operator=(const TransportCableModel &) = delete;
  TransportCableModel &operator=(TransportCableModel &&) = delete;

  void clear();
  int validate();
  std::vector<Pile> m_pilesUp = {};
  std::vector<Pile> m_pilesDown = {};
};
}  // namespace siconos::fem::cable
