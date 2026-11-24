// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "CableTools.h"

#include <fstream>

nlohmann::json siconos::fem::cable::tools::load_json_file(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file '" + filename + "'");
  }

  nlohmann::json j;
  file >> j;
  return j;
}
