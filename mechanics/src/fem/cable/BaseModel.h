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
  BaseModel(const BaseModel &) = default;
  BaseModel &operator=(const BaseModel &) = default;
  BaseModel(BaseModel &&) = default;
  BaseModel &operator=(BaseModel &&) = default;
  virtual ~BaseModel() noexcept = default;

  void from_json(const json &j, const std::string &a_header);

 protected:
  virtual void from_json(const json &j) = 0;
};
}  // namespace siconos::fem::cable
#endif
