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

/*! \file SiconosVisitor.hpp
  \brief A general visitor interface for siconos objects.
*/

#ifndef SiconosVisitor_hpp
#define SiconosVisitor_hpp

/** A visitor pattern.

   User has to instantiate a derivation of SiconosVisitor class :

   struct myvisitor : public SiconosVisitor
   {
     void visit(const LagrangianDS& ds)
     {
       ...
     }
   }

   with some wanted visit() functions.

   Then the visitor may be used as :

   A_visitable_Siconos_Object->accept(Siconos::Visitor myvisitor)

   SiconosVisitor also define a type visitor object under the
   namespace Type:: and some functions to access type of visitables
   classes:

   Type::value(c) : the type of the visitable object c as an enum.

   Type::name(c)  : the name of the Type::value as a std::string

*/

#include <memory>

#include "SiconosException.hpp"

namespace siconos::algebra {
class SimpleMatrix;
class BlockMatrix;
class SiconosVector;
class BlockVector;
}  // namespace siconos::algebra

// We have to declare all the classes that might be visited.
namespace siconos::modeling {
// NS laws
class EqualityConditionNSL;
class ComplementarityConditionNSL;
class MixedComplementarityConditionNSL;
class MultipleImpactNSL;
class NewtonImpactFrictionNSL;
class NewtonImpactNSL;
class NewtonImpactRollingFrictionNSL;
class NormalConeNSL;
class RelayNSL;
class FremondImpactFrictionNSL;

// DS
class FirstOrderNonLinearDS;
class NewtonEulerDS;
class LagrangianDS;

// Relations
class LagrangianR;
class NewtonEulerR;
class FirstOrderR;
}  // namespace siconos::modeling

namespace siconos::internal {

struct SiconosVisitor {
  SiconosVisitor() = default;
  SiconosVisitor(const SiconosVisitor &) = delete;
  SiconosVisitor(SiconosVisitor &&) = delete;
  SiconosVisitor &operator=(SiconosVisitor &&) = delete;
  SiconosVisitor &operator=(const SiconosVisitor &) = delete;

  // declaration of visitors used in Interactions to set level values.
  // Note FP: use variant for nslaws ?
  // Pros: less virtual functions (no more class hierarchy in nslaws), can use std::visit
  // Cons: must use std::visit. Won't change anything at user level, but it will be a bit
  // more difficult to call functions at low level.
  virtual void visit(const siconos::modeling::ComplementarityConditionNSL &nslaw) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };
  virtual void visit(const siconos::modeling::EqualityConditionNSL &nslaw) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };
  virtual void visit(const siconos::modeling::MixedComplementarityConditionNSL &nslaw) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };
  virtual void visit(const siconos::modeling::MultipleImpactNSL &nslaw) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };
  virtual void visit(const siconos::modeling::NewtonImpactFrictionNSL &nslaw) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };
  virtual void visit(const siconos::modeling::NewtonImpactNSL &nslaw) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };
  virtual void visit(const siconos::modeling::NewtonImpactRollingFrictionNSL &nslaw) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };
  virtual void visit(const siconos::modeling::NormalConeNSL &nslaw) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };
  virtual void visit(const siconos::modeling::RelayNSL &nslaw) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };
  virtual void visit(const siconos::modeling::FremondImpactFrictionNSL &nslaw) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };

  // Dynamical Systems
  virtual void visit(const siconos::modeling::LagrangianDS &ds) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };

  virtual void visit(const siconos::modeling::FirstOrderNonLinearDS &ds) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };

  virtual void visit(const siconos::modeling::NewtonEulerDS &ds) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };

  // Relations
  virtual void visit(const siconos::modeling::LagrangianR &) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };

  virtual void visit(const siconos::modeling::FirstOrderR &) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };

  virtual void visit(const siconos::modeling::NewtonEulerR &) const {
    THROW_EXCEPTION("you must define a visit function in a derived class of SiconosVisitor");
  };

  virtual ~SiconosVisitor() noexcept = default;
};

}  // namespace siconos::internal

#endif /* SiconosVisitor_hpp */
