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

// /** Print the contents of any sequence - From Thinking in C++, vol 2 p365.
// \param first, any iterator, beginning of the sequence
// \param last, any iterator, end of the sequence
// \param char*, optional message on top of output, default ""
// \param char*, separator between sequence elements, default new line
// \param ostream, output destination, default std::cout
// */
// template <typename Iter>
// void print(Iter first, Iter last, const char* nm = "", const char* sep = "\n",
//            std::ostream& os = std::cout)
// {
//   if (nm != nullptr && *nm != '\0') os << nm << ": " << sep;
//   typedef typename std::iterator_traits<Iter>::value_type T;
//   std::copy(first, last, std::ostream_iterator<T>(os, sep));
//   os << std::endl;
// }

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

}  // namespace siconos::tools

#endif
