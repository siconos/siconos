/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "LagrangianSparseR.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianSparseDS.hpp"
#include "Tools.hpp"

void siconos::modeling::LagrangianSparseR::display() const {
  std::cout << "=====> Relation of type " << tools::enum_to_string(_relationType)
            << " and subtype " << tools::enum_to_string(_subType) << "\n";
  if (hasJacobianhOver_q_) {
    std::cout << "- jacobian of h over q\n";
    useJacobianhOver_q([&](const auto& M) { siconos::algebra::print(M); });
    std::cout << "\n";
  }
}

void siconos::modeling::LagrangianSparseR::allocate_read_dynamical_systems_var_vectors(
    std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars,
    modeling::DynamicalSystem& ds1, modeling::DynamicalSystem& ds2) const {
  ds_vars.resize(tools::enum_to_index(LagrangianSparseR::ds_var::size));

  // Default ds_vars
  ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q0)] =
      std::make_shared<siconos::algebra::BlockVector>();
  ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q1)] =
      std::make_shared<siconos::algebra::BlockVector>();

  bool has_two_ds = &ds1 != &ds2;

  if (auto* lag1 = dynamic_cast<LagrangianDS*>(&ds1)) {
    // Put q, velocity of each DS into a block. (Pointers links, no copy!!)
    ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q0)]->insertPtr(lag1->q());
    ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q1)]->insertPtr(lag1->velocity());

    if (lag1->acceleration()) {
      if (!ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)])
        ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)] =
            std::make_shared<siconos::algebra::BlockVector>();  // acceleration

      ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)]->insertPtr(
          lag1->acceleration());
    }

    for (unsigned int k = 0; k < 3; k++) {
      if (lag1->p(k)) {
        if (!ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k])
          ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k] =
              std::make_shared<siconos::algebra::BlockVector>();
        ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k]->insertPtr(
            lag1->p(k));
      }
    }
  } else if (auto* lag1 = dynamic_cast<LagrangianSparseDS*>(&ds1)) {
    // Put q, velocity of each DS into a block. (Pointers links, no copy!!)
    ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q0)]->insertPtr(lag1->q());
    ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q1)]->insertPtr(lag1->velocity());

    if (lag1->acceleration()) {
      if (!ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)])
        ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)] =
            std::make_shared<siconos::algebra::BlockVector>();  // acceleration

      ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)]->insertPtr(
          lag1->acceleration());
    }

    for (unsigned int k = 0; k < 3; k++) {
      if (lag1->p(k)) {
        if (!ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k])
          ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k] =
              std::make_shared<siconos::algebra::BlockVector>();
        ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k]->insertPtr(
            lag1->p(k));
      }
    }
  }

  if (has_two_ds) {
    if (auto* lag2 = dynamic_cast<LagrangianDS*>(&ds2)) {
      // Put q, velocity of each DS into a block. (Pointers links, no copy!!)
      ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q0)]->insertPtr(lag2->q());
      ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q1)]->insertPtr(
          lag2->velocity());

      if (lag2->acceleration()) {
        if (!ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)])
          ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)] =
              std::make_shared<siconos::algebra::BlockVector>();  // acceleration

        ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)]->insertPtr(
            lag2->acceleration());
      }

      for (unsigned int k = 0; k < 3; k++) {
        if (lag2->p(k)) {
          if (!ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k])
            ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k] =
                std::make_shared<siconos::algebra::BlockVector>();
          ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k]->insertPtr(
              lag2->p(k));
        }
      }
    } else if (auto* lag2 = dynamic_cast<LagrangianSparseDS*>(&ds2)) {
      // Put q, velocity of each DS into a block. (Pointers links, no copy!!)
      ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q0)]->insertPtr(lag2->q());
      ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q1)]->insertPtr(
          lag2->velocity());

      if (lag2->acceleration()) {
        if (!ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)])
          ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)] =
              std::make_shared<siconos::algebra::BlockVector>();  // acceleration

        ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::q2)]->insertPtr(
            lag2->acceleration());
      }

      for (unsigned int k = 0; k < 3; k++) {
        if (lag2->p(k)) {
          if (!ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k])
            ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k] =
                std::make_shared<siconos::algebra::BlockVector>();
          ds_vars[tools::enum_to_index(LagrangianSparseR::ds_var::p0) + k]->insertPtr(
              lag2->p(k));
        }
      }
    }
  }
}

