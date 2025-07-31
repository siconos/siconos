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
#include "TransportCableResult.h"

#include <fstream>
#include <memory>

#include "PulleyWrapping.h"
#include "Rope.h"  // IWYU pragma: keep
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"  // IWYU pragma: keep

siconos::fem::cable::TransportCableResult::TransportCableResult() {
  ropes_down.set_Down(true);
}

void siconos::fem::cable::TransportCableResult::prepareSupport() {
  supports.clear();
  topPulleyId = -1;
  downPulleyId = -1;
  ropes_up.prepareSupport(supports, downPulleyId);
  ropes_down.prepareSupport(supports, topPulleyId);

  // Bottom pylon/pulley
  supports[downPulleyId]->prepare(ropes_down.getFirstPylon(), ropes_up.getFirstPylon(),
                                  ropes_down.initialTension());

  // Top pylon/pulley
  supports[topPulleyId]->prepare(ropes_up.getLastPylon(), ropes_down.getLastPylon(),
                                 ropes_up.getTensionAtLastNode());
}

void siconos::fem::cable::TransportCableResult::prepareIneqConstraint(int nb_node) {
  contacts.clear();
  contacts.resize(nb_node, 0);
  gVector.resize(nb_node);
  gVector.setConstant(1.);
  // jacobian_g_Over_q =
  //     std::make_shared<siconos::algebra::SiconosDenseMatrix>(nb_node, 3 * nb_node);
  jacobian_g_Over_q.resize(3 * nb_node, nb_node);

  T = std::make_shared<siconos::algebra::SiconosDenseMatrix>(nb_node, 3 * nb_node);
  jacobian_g_Over_q.setZero();
  T->setZero();
}

int siconos::fem::cable::TransportCableResult::exportTC(const std::string &a_fileName,
                                                        nlohmann::ordered_json &a_output,
                                                        const std::string &a_option) {
  int res = EXIT_SUCCESS;
  try {
    res = to_json(a_output, a_option);
  } catch (const nlohmann::json::exception &ex) {
    // "Error exporting model " << ex.what());
    throw ex;
  }
  if (a_fileName != "") {
    res = EXIT_FAILURE;
    std::ofstream outputFile(a_fileName);
    if (outputFile.is_open()) {
      outputFile << a_output.dump(2);
      res = EXIT_SUCCESS;
    }
  }
  return res;
}

int siconos::fem::cable::TransportCableResult::to_json(nlohmann::ordered_json &j,
                                                       const std::string &a_option) {
  if (a_option == "fem") {
    j["g"] = gVector;
  } else if (a_option == "ropeway") {
    j["ropes_up"] = ropes_up;
    j["ropes_down"] = ropes_down;

    // ropes_up.to_json(j["ropes_up"]);
    // ropes_down.to_json(j["ropes_down"]);
  } else {
    j["ropes_up"] = ropes_up;
    j["ropes_down"] = ropes_down;

    // ropes_up.to_json(j["ropes_up"]);
    // ropes_down.to_json(j["ropes_down"]);
    j["supports"] = nlohmann::ordered_json::array();
    j["pulleys"] = nlohmann::ordered_json::array();
    for (auto s : supports) {
      if (auto pulley = std::dynamic_pointer_cast<PulleyWrapping>(s))
        j["pulleys"].push_back(*pulley);
      else
        j["supports"].push_back(*s);
    }
    if (q.size() > 0) {
      j["q"] = q;
      j["R"] = R;
      j["tension"] = TS;
      j["punct"] = weightVector;
      j["cable_length"] = totalLength;
      j["contacts"] = contacts;
      j["g"] = gVector;
    }
  }
  return EXIT_SUCCESS;
}
