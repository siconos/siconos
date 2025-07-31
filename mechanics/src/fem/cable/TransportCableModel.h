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
 *  The model can be:
 *  - build empty and so must be read from an input file using from_json or from_file
 *  - build directlty from a file or json input
 *  Read-only use afterwards.
 *
 *  An example of the format of this input file is
 *  available in mechanics/src/fem/cable/test/data/model.origin.json
 *
 */
class TransportCableModel {
 private:
  MechanicalProperties mechanicalProperties_;

  /** The set of vehicles */
  Carriers m_carriers;

  /** Pylons corresponding to the 'up' ropeway, including end-line stations */
  std::vector<Pylon> list_of_pylons_up_ = {};

  /** Pylons corresponding to the 'down' ropeway, including end-line stations */
  std::vector<Pylon> list_of_pylons_down_ = {};

  // Rule of five
  TransportCableModel(const TransportCableModel &) = delete;
  TransportCableModel(TransportCableModel &&) = delete;
  TransportCableModel &operator=(const TransportCableModel &) = delete;
  TransportCableModel &operator=(TransportCableModel &&) = delete;

 public:
  TransportCableModel() = delete;

  virtual ~TransportCableModel() noexcept = default;

  /** Read data from a json input
      \param input json object
  */
  explicit TransportCableModel(const nlohmann::json &input);

  /** \return true if the model has been loaded and validated */
  bool isLoaded() const;

  /** \return the mechanical properties of the cable */
  inline const MechanicalProperties &mechanicalProperties() const {
    return mechanicalProperties_;
  }

  /** \return the set of vehicles */
  inline const Carriers &get_carriers() const { return m_carriers; };

  /** \return the list of pylons corresponding to the 'up' ropeway */
  inline const std::vector<Pylon> &get_pylons_up() const { return list_of_pylons_up_; };

  /** \return the list of pylons corresponding to the 'down' ropeway */
  inline const std::vector<Pylon> &get_pylons_down() const { return list_of_pylons_down_; };
};

}  // namespace siconos::fem::cable

// Serialization
namespace nlohmann {
template <>
struct adl_serializer<siconos::fem::cable::TransportCableModel> {
  static siconos::fem::cable::TransportCableModel from_json(const json &j) {
    return siconos::fem::cable::TransportCableModel(j);
  }
  // static void to_json(json &j, const siconos::fem::cable::TransportCableModel &input) {
  //   j = json{};
  // }
};
}  // namespace nlohmann