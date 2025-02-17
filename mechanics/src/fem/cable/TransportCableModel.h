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
  Class TransportCableModel : full description of the setup (cable, piles, carriers)
  before discretization
*/
#pragma once

#include <vector>

#include "Carriers.h"
#include "MechanicalProperties.h"
#include "Pylon.h"

namespace siconos::fem::cable {

/** A class to handle the full description of the setup
 *
 *  - mechanical properties
 *  - the pylons (standard + top and down stations)
 *  - carriers
 *
 *  Build from json and read-only use afterwards.
 *
 */
class TransportCableModel {
 private:
  MechanicalProperties m_cable;

  /** The set of vehicles */
  Carriers m_carriers;

  /** End station (top) */
  Pylon m_stationUp;

  /** Down station */
  Pylon m_stationDown;

  /** Standard pylons list (without top and down station) */
  std::vector<Pylon> m_piles = {};

  /** Pylons corresponding to the 'up' ropeway, including end-line stations */
  std::vector<Pylon> m_pilesUp = {};

  /** Pylons corresponding to the 'down' ropeway, including end-line stations */
  std::vector<Pylon> m_pilesDown = {};

  // Rule of five
  TransportCableModel(const TransportCableModel &) = delete;
  TransportCableModel(TransportCableModel &&) = delete;
  TransportCableModel &operator=(const TransportCableModel &) = delete;
  TransportCableModel &operator=(TransportCableModel &&) = delete;

  /** Reset the model */
  void clear();

  /** Internal validation of the model (fill pilesUp and pilesDown in) */
  int validate();

 public:
  TransportCableModel() = default;

  /** Read data from a json file
      \param a_filename input file name (json)
  */
  TransportCableModel(const std::string &a_filename);

  virtual ~TransportCableModel() noexcept = default;

  /** Build the model from a file (json)
   * \param a_filename json file name
   */
  int from_file(const std::string &a_fileName);

  /** Build a model from a json object
   *  \param j json object
   */
  int from_json(const nlohmann::json &j);

  /** \return true if the model has been loaded and validated */
  bool isLoaded();

  /** \return the mechanical properties of the cable */
  const MechanicalProperties &mechanicalProperties() const;

  /** \return the set of vehicles */
  const Carriers &get_carriers() const;

  /** \return the list of pylons corresponding to the 'up' ropeway */
  const std::vector<Pylon> &get_pylons_up() const;

  /** \return the list of pylons corresponding to the 'down' ropeway */
  const std::vector<Pylon> &get_pylons_down() const;
};
}  // namespace siconos::fem::cable
