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

/*! \file Point.h
  Geometrical point

*/
#pragma once

#include "BaseModel.h"

namespace siconos::fem::cable {
class Point : public BaseModel {
 public:
  double x{0.};
  double y{0.};
  double z{0.};

  Point() = default;
  Point(const Point &) = default;
  Point(Point &&) = default;
  Point &operator=(const Point &) = default;
  Point &operator=(Point &&) = default;

  Point(double ax, double ay, double az) : x{ax}, y{ay}, z{az} {};
  virtual ~Point() noexcept = default;

  virtual void from_json(const json &j);

  double norm() { return sqrt(x * x + y * y + z * z); }

  void diff(const Point &p1, const Point &p2) {
    x = p1.x - p2.x;
    y = p1.y - p2.y;
    z = p1.z - p2.z;
  }
  void add(const Point &p1, const Point &p2) {
    x = p1.x + p2.x;
    y = p1.y + p2.y;
    z = p1.z + p2.z;
  }
  void mult(double a_scalar) {
    x *= a_scalar;
    y *= a_scalar;
    z *= a_scalar;
  }
  void opposite() {
    x = -x;
    y = -y;
    z = -z;
  }
  void zero() {
    x = 0;
    y = 0;
    z = 0;
  }
  double dot(const Point &p) { return x * p.x + y * p.y + z * p.z; }
};

bool operator>(const Point &p1, const Point &p2);
bool operator<(const Point &p1, const Point &p2);
void to_json(ojson &j, const Point &p);

}  // namespace siconos::fem::cable
