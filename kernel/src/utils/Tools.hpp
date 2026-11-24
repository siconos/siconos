// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file Tools.hpp
  Various useful functions and typedef.
*/

#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <limits>
#include <ranges>
#include <sstream>

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
