/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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
#include <nlohmann/json.hpp>
#include <vector>

namespace siconos::fem::cable {

class Point;
}

namespace siconos::fem::cable::tools {

// Note FP: use a generic name. We may decide later to replace json by hdf5 or something else
/** Type of input used to read data */
using InputNodeType = nlohmann::json;

/** Type of output used to write data */
using OutpuNodeType = nlohmann::ordered_json;

template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T, std::integral S>
std::vector<T> linspace(T start, T end, S num) {
  std::vector<T> result(num);

  if (num == 0) {
    return result;
  }

  T step = (end - start) / (num - 1);

  for (S i = 0; i < num; ++i) {
    result[i] = (start + step * i);
  }
  result[num - 1] = end;  // I want to ensure that start and end
                          // are exactly the same as the input
  return result;
}

/**
 * Retrieve a parameter from a json object
 *
 * This template function extracts a parameter of type `T` from a json object.
 * If the json object does not contain the given field name or is null,
 * the provided default value is returned.
 *
 * \tparam T Type of the parameter to extract (must be deserializable from JSON).
 * \param a_arg The JSON object to read from.
 * \param a_name The name of the field to look for in the JSON object.
 * \param a_default The default value to return if the field is not found or JSON is null.
 * \return T The extracted value from the JSON or the default if not present.
 *
 */
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

/** Read a json file

  \param filename name of the input file
  \return a json object
*/
nlohmann::json load_json_file(const std::string &filename);

}  // namespace siconos::fem::cable::tools

#endif
