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
#include "TransportCableModel.h"

#include <algorithm>
#include <cstdlib>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

siconos::fem::cable::TransportCableModel::TransportCableModel(const nlohmann::json &input)
    : m_carriers(input.at("carriers")) {
  // -- Read cable (mechanical properties) and carriers (vehicle positions) --
  if (!input.contains("mechanicalProperties")) throw std::runtime_error("Missing 'cable' key in json input.");

  if (!input.contains("carriers"))
    throw std::runtime_error("Missing 'carriers' key in json input.");

  if (!input.contains("stationUp"))
    throw std::runtime_error("Missing 'stationUp' key in json input.");
  if (!input.contains("stationDown"))
    throw std::runtime_error("Missing 'stationDown' key in json input.");

  if (!input.contains("piles")) throw std::runtime_error("Missing 'piles' key in json input.");

  // cable mechanical properties
  mechanicalProperties_ = input["mechanicalProperties"];

  list_of_pylons_up_.clear();
  list_of_pylons_down_.clear();

  // Read stations positions
  Pylon stationUp(input["stationUp"], true);
  Pylon stationDown(input["stationDown"], true);
  assert(stationDown < stationUp);

  // -- Read pylons positions (excluding station up and down) --
  std::vector<Pylon> list_of_pylons = {};
  const json &jpiles = input["piles"];
  assert(jpiles.is_array());
  list_of_pylons.reserve(jpiles.size());
  for (const auto &jp : jpiles) {
    list_of_pylons.emplace_back(jp, false);
  }
  assert(list_of_pylons.size());
  // piles x-coords must be in ascending order
  std::sort(list_of_pylons.begin(), list_of_pylons.end());

  // Ensure stations are at correct positions
  assert((stationDown < list_of_pylons.front() && list_of_pylons.back() < stationUp));

  // Build the list corresponding to "up" pylons (including stations) from list_of_pylons
  list_of_pylons_up_.emplace_back(stationDown);
  for (auto &p : list_of_pylons) {
    list_of_pylons_up_.push_back(p);  // copy
  }
  list_of_pylons_up_.emplace_back(stationUp);

  // Now the list of "down" pylons (no stations) with a shift corresponding to the distance
  // between the ropes
  for (auto &p : list_of_pylons_up_) {
    list_of_pylons_down_.push_back(p);  // copy
    list_of_pylons_down_.back().shift_y();
  }

  if (!isLoaded())
    throw std::runtime_error(
        "Something went wrong during model loading: there are no pylons (empty lists).");
}

bool siconos::fem::cable::TransportCableModel::isLoaded() const {
  return (list_of_pylons_up_.size() != 0 && list_of_pylons_down_.size() != 0);
}