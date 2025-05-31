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

#ifndef CABLETOOLS
#define CABLETOOLS
#include <memory>
#include <nlohmann/json.hpp>
#include <vector>

#include "SiconosVector.hpp"

namespace siconos::fem::cable {

class Point;
}

namespace siconos::fem::cable::tools {

void pointsToSiconosVector(const std::vector<siconos::fem::cable::Point> vecin,
                           std::shared_ptr<siconos::algebra::SiconosVector> vecout);

template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in) {
  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) {
    return linspaced;
  }
  if (num == 1) {
    linspaced.push_back(start);
    return linspaced;
  }

  double delta = (end - start) / (num - 1);

  for (int i = 0; i < num - 1; ++i) {
    linspaced.push_back(start + delta * i);
  }
  linspaced.push_back(end);  // I want to ensure that start and end
                             // are exactly the same as the input
  return linspaced;
}

template <typename T>
T getParam(const nlohmann::json &a_arg, const std::string &a_name, T a_default) {
  T vRet = a_default;
  if (!a_arg.is_null()) {
    if (a_arg.contains(a_name)) {
      a_arg.at(a_name).get_to(vRet);
    }
  }
  return vRet;
}

}  // namespace siconos::fem::cable::tools

#endif
