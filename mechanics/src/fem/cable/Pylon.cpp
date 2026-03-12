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

#include "Pylon.h"

#include "SiconosMatrix.hpp"

void siconos::fem::cable::Pylon::from_json(const nlohmann::json &j) {
  coordinates_ << j.at("x"), j.at("y"), j.at("z");
  j.at("R").get_to(radius_);
  distanceToUpRope_ = j.value("dUp", 0.);  // 0 default for station pylons
  distanceToDownRope_ = j.value("dDown", 0.);
  height_ = j.value("h", 0.);  // h is not always provided
}

void siconos::fem::cable::Pylon::shift_y() {
  if (isAStation_) {
    coordinates_(1) += 2.0 * radius_;
  } else {
    coordinates_(1) += distanceToUpRope_ + distanceToDownRope_;
  }
}

void siconos::fem::cable::Pylon::display() const {
  std::cout << "--- Pylon: \n";
  if (isAStation_)
    std::cout << "- station type\n";
  else
    std::cout << "- standard type\n";
  std::cout << "- position: ";
  siconos::algebra::print(coordinates_.transpose());
  std::cout << "- radius: " << radius_;
  std::cout << " - heigth: " << height_;
  std::cout << " - distance to up and down ropes: " << distanceToUpRope_ << ","
            << distanceToDownRope_ << "\n --------------\n\n";
}

bool siconos::fem::cable::operator<(const siconos::fem::cable::Pylon &p1,
                                    const siconos::fem::cable::Pylon &p2) {
  return (p1.coords()(0) < p2.coords()(0));
}
