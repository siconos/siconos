// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file BaseModel.h
  \brief
*/

#ifndef BASEMODEL_H
#define BASEMODEL_H

#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;
using ojson = nlohmann::ordered_json;

namespace siconos::fem::cable {

class BaseModel {
 public:
  // Rule of 5
  BaseModel() = default;
  BaseModel(const BaseModel&) = default;
  BaseModel& operator=(const BaseModel&) = default;
  BaseModel(BaseModel&&) = default;
  BaseModel& operator=(BaseModel&&) = default;
  virtual ~BaseModel() noexcept = default;

  void from_json(const json& j, const std::string& a_header);

 protected:
  virtual void from_json(const json& j) = 0;
};
}  // namespace siconos::fem::cable
#endif
