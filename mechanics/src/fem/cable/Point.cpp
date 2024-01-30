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
#include "Point.h"

#include <iostream>

void siconos::fem::cable::Point::from_json(const json &j) {
  j.at("x").get_to(x);
  j.at("y").get_to(y);
  j.at("z").get_to(z);
}

bool siconos::fem::cable::operator>(const Point &p1, const Point &p2) {
  return (p1.x >= p2.x);
}

bool siconos::fem::cable::operator<(const Point &p1, const Point &p2) { return (p1.x < p2.x); }

void siconos::fem::cable::to_json(ojson &j, const Point &p) {
  j.push_back(p.x);
  j.push_back(p.y);
  j.push_back(p.z);
}
