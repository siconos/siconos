/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "BaseModel.h"

void siconos::fem::cable::BaseModel::from_json(const json& j, const std::string& a_header) {
  if (j.contains(a_header)) {
    from_json(j[a_header]);
  }
}
