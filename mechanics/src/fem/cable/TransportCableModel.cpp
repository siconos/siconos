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
#include <fstream>
#include <iostream>
#include <map>
#include <nlohmann/json.hpp>
#include <string>

#include "Cable.h"

using json = nlohmann::json;

siconos::fem::cable::TransportCableModel::TransportCableModel(const std::string &a_filename) {
  from_file(a_filename);
}

int siconos::fem::cable::TransportCableModel::from_file(const std::string &a_fileName) {
  int res = EXIT_FAILURE;
  json input;
  std::ifstream file(a_fileName);
  if (file.is_open()) {
    file >> input;
    res = from_json(input);
  }
  return res;
}

int siconos::fem::cable::TransportCableModel::from_json(const json &j) {
  int res = EXIT_SUCCESS;
  std::map<std::string, BaseModel &> vElems = {{"cable", m_cable},
                                               {"carriers", m_carriers},
                                               {"stationUp", m_stationUp},
                                               {"stationDown", m_stationDown}};
  clear();
  try {
    for (auto &vElem : vElems) {
      vElem.second.from_json(j, vElem.first);
    }
    if (j.contains("piles")) {
      const json &jpiles = j["piles"];
      if (jpiles.is_array()) {
        m_piles.reserve(jpiles.size());
        for (const auto &jp : jpiles) {
          m_piles.push_back(Pile());
          m_piles.back().from_json(jp);
        }
      }
    }
    res = validate();
  } catch (const json::exception &) {
    res = EXIT_FAILURE;
  }

  return res;
}

bool siconos::fem::cable::TransportCableModel::isLoaded() {
  return (m_pilesUp.size() != 0 && m_pilesDown.size() != 0);
}

const siconos::fem::cable::Cable &siconos::fem::cable::TransportCableModel::get_cable() const {
  return m_cable;
}

const siconos::fem::cable::Carriers &siconos::fem::cable::TransportCableModel::get_carriers()
    const {
  return m_carriers;
}

const std::vector<siconos::fem::cable::Pile> &
siconos::fem::cable::TransportCableModel::get_piles1() const {
  return m_pilesUp;
}

const std::vector<siconos::fem::cable::Pile> &
siconos::fem::cable::TransportCableModel::get_piles2() const {
  return m_pilesDown;
}

void siconos::fem::cable::TransportCableModel::clear() {
  // raz du modèle
  m_piles.clear();
  m_pilesUp.clear();
  m_pilesDown.clear();
}

int siconos::fem::cable::TransportCableModel::validate() {
  // Validation du modèle
  // Création des 2 lignes à partir de la définition des poteaux
  if (m_stationUp > m_stationDown) {
    bool vOk = true;
    if (m_piles.size()) {
      // les x des poteaux doivent être croissant
      std::sort(m_piles.begin(), m_piles.end());

      vOk = (m_piles.front() > m_stationDown && m_piles.back() < m_stationUp);
    }
    if (vOk) {
      m_pilesUp.push_back(Pile(m_stationDown, true));
      for (auto &p : m_piles) {
        m_pilesUp.push_back(p);  // Copy
      }
      m_pilesUp.push_back(Pile(m_stationUp, true));

      for (auto &p : m_pilesUp) {
        m_pilesDown.push_back(p);  // Copy
        p.transform(true);         // useless ???
        m_pilesDown.back().transform(false);
      }
      return EXIT_SUCCESS;
    }
  }

  return EXIT_FAILURE;
}
