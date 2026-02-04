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

/*! \file Tools.hpp
  Various useful functions and typedef.
*/

#ifndef TOOLS_H
#define TOOLS_H

// #include <algorithm>
#include <iostream>
// #include <iterator>
#include <ranges>
#include <sstream>
// #include <string>

namespace siconos::tools {

constexpr auto PBSTR = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
constexpr auto PBWIDTH = 60;

/** A function to convert any type to std::string*/
template <class T>
std::string toString(const T& obj) {
  static std::ostringstream o;
  o.str("");
  o << obj;
  return o.str();
}

/** Print to screen the content of any stl-like container.
    Warning: c++20 required
    \param description anything like a string, a number ... to describe the output
    \param container the object which content must be printed
*/
void print(auto const description, std::ranges::range auto const& container) {
  for (std::cout << description; auto const& e : container) std::cout << e << ' ';
  std::cout << '\n';
}

void progressBar(double percentage);

template <typename T>
std::string enum_to_string(T value) {
  return std::to_string(static_cast<typename std::underlying_type<T>::type>(value));
}

// Helper to access to the underlying index of an enum class
template <typename Enum>
constexpr auto enum_to_index(Enum e) noexcept {
  return static_cast<std::underlying_type_t<Enum>>(e);
}

}  // namespace siconos::tools

#endif
