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

/*! \file TypeName.hpp
  \brief get a string name from visitable classes
*/

#ifndef TypeName_hpp
#define TypeName_hpp

#include <cassert>
#include <string>
#include <variant>

namespace siconos::modeling {

class FirstOrderNonLinearDS;
class NewtonEulerDS;
class LagrangianDS;
class FirstOrderLinearDS;
class FirstOrderLinearTIDS;
class LagrangianLinearTIDS;
class LagrangianLinearDiagonalDS;
class FremondImpactFrictionNSL;
class RelayNSL;
class NormalConeNSL;
class NewtonImpactRollingFrictionNSL;
class NewtonImpactNSL;
class NewtonImpactFrictionNSL;
class MultipleImpactNSL;
class MixedComplementarityConditionNSL;
class EqualityConditionNSL;
class ComplementarityConditionNSL;

enum class Type {
  FirstOrderNonLinearDS,
  NewtonEulerDS,
  LagrangianDS,
  FirstOrderLinearDS,
  FirstOrderLinearTIDS,
  LagrangianLinearTIDS,
  LagrangianLinearDiagonalDS,
  // NSLaws
  RelayNSL,
  NormalConeNSL,
  NewtonImpactRollingFrictionNSL,
  NewtonImpactNSL,
  FremondImpactFrictionNSL,
  NewtonImpactFrictionNSL,
  MultipleImpactNSL,
  MixedComplementarityConditionNSL,
  EqualityConditionNSL,
  ComplementarityConditionNSL,
};

}  // namespace siconos::modeling

namespace siconos::types {

struct FindType {
  // auto operator()(const siconos::experimental::B& obj) const {return experimental::Type::B;
  // } auto operator()(const siconos::experimental::C& obj) const {return
  // experimental::Type::C; } auto operator()(std::shared_ptr<siconos::experimental::B> obj)
  // const {return experimental::Type::B; } auto
  // operator()(std::shared_ptr<siconos::experimental::C> obj) const {return
  // experimental::Type::C; }

  // DS
  auto visit(const siconos::modeling::FirstOrderNonLinearDS&) const {
    return siconos::modeling::Type::FirstOrderNonLinearDS;
  };
  auto visit(const siconos::modeling::NewtonEulerDS&) const {
    return siconos::modeling::Type::NewtonEulerDS;
  };
  auto visit(const siconos::modeling::LagrangianDS&) const {
    return siconos::modeling::Type::LagrangianDS;
  };
  auto visit(const siconos::modeling::LagrangianLinearTIDS&) const {
    return siconos::modeling::Type::LagrangianLinearTIDS;
  };
  auto visit(const siconos::modeling::LagrangianLinearDiagonalDS&) const {
    return siconos::modeling::Type::LagrangianLinearDiagonalDS;
  };

  // NSlaws
  auto visit(const siconos::modeling::RelayNSL&) const {
    return siconos::modeling::Type::RelayNSL;
  };
  auto visit(const siconos::modeling::NormalConeNSL&) const {
    return siconos::modeling::Type::NormalConeNSL;
  };
  auto visit(const siconos::modeling::NewtonImpactRollingFrictionNSL&) const {
    return siconos::modeling::Type::NewtonImpactRollingFrictionNSL;
  };
  auto visit(const siconos::modeling::NewtonImpactNSL&) const {
    return siconos::modeling::Type::NewtonImpactNSL;
  };
  auto visit(const siconos::modeling::NewtonImpactFrictionNSL&) const {
    return siconos::modeling::Type::NewtonImpactFrictionNSL;
  };
  auto visit(const siconos::modeling::MultipleImpactNSL&) const {
    return siconos::modeling::Type::MultipleImpactNSL;
  };
  auto visit(const siconos::modeling::MixedComplementarityConditionNSL&) const {
    return siconos::modeling::Type::MixedComplementarityConditionNSL;
  };
  auto visit(const siconos::modeling::EqualityConditionNSL&) const {
    return siconos::modeling::Type::EqualityConditionNSL;
  };
  auto visit(const siconos::modeling::FremondImpactFrictionNSL&) const {
    return siconos::modeling::Type::MixedComplementarityConditionNSL;
  };
  auto visit(const siconos::modeling::ComplementarityConditionNSL&) const {
    return siconos::modeling::Type::ComplementarityConditionNSL;
  };
};

// auto type_value = []( auto& obj) { return obj.name(); };

static FindType find;

template <typename C>
inline auto type_value(const C& c) {
  return c.acceptType(find);
}

template <typename T>
constexpr auto str(const T& X) {
  switch (X) {
    case T::FirstOrderNonLinearDS:
      return "siconos::modeling::FirstOrderNonLinearDS";
      break;
    case T::NewtonEulerDS:
      return "siconos::modeling::NewtonEulerDS";
      break;
    case T::LagrangianDS:
      return "siconos::modeling::LagrangianDS";
      break;
    case T::LagrangianLinearTIDS:
      return "siconos::modeling::LagrangianLinearTIDS";
      break;
    case T::LagrangianLinearDiagonalDS:
      return "siconos::modeling::LagrangianLinearDiagonalDS";
      break;
    case T::RelayNSL:
      return "siconos::modeling::RelayNSL";
      break;
    case T::FremondImpactFrictionNSL:
      return "siconos::modeling::FremondImpactFrictionNSL";
      break;
    case T::NormalConeNSL:
      return "siconos::modeling::NormalConeNSL";
      break;
    case T::NewtonImpactRollingFrictionNSL:
      return "siconos::modeling::NewtonImpactRollingFrictionNSL";
      break;
    case T::NewtonImpactNSL:
      return "siconos::modeling::NewtonImpactNSL";
      break;
    case T::NewtonImpactFrictionNSL:
      return "siconos::modeling::NewtonImpactFrictionNSL";
      break;
    case T::MultipleImpactNSL:
      return "siconos::modeling::MultipleImpactNSL";
      break;
    case T::MixedComplementarityConditionNSL:
      return "siconos::modeling::MixedComplementarityConditionNSL";
      break;
    case T::EqualityConditionNSL:
      return "siconos::modeling::EqualityConditionNSL";
      break;
    case T::ComplementarityConditionNSL:
      return "siconos::modeling::ComplementarityConditionNSL";
      break;
    default:
      assert(0);
      return "unknown";
  }
}

template <class C>
std::string type_name(const C& c) {
  return str(type_value(c));
}

}  // namespace siconos::types

#endif
