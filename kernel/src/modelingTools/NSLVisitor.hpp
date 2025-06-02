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

/*! \file NSLVisitor.hpp
  Abstract class - General interface for visitors of  nonsmooth laws.
*/

#ifndef NSLVisitor_hpp
#define NSLVisitor_hpp

/** A visitor pattern.

   User has to instantiate a derivation of NSLVisitor class :

   struct myvisitor : public NSLVisitor
   {
     void visit(const LagrangianDS& ds)
     {
       ...
     }
   }

   with some wanted visit() functions.

   Then the visitor may be used as :

   A_visitable_Siconos_Object->accept(siconos::Visitor myvisitor)

   NSLVisitor also define a type visitor object under the
   namespace Type:: and some functions to access type of visitables
   classes:

   Type::value(c) : the type of the visitable object c as an enum.

   Type::name(c)  : the name of the Type::value as a std::string

*/

#include <memory>

#include "SiconosConfig.h"  // For HAVE_SICONOS ...
#include "SiconosException.hpp"

// We have to declare all the classes that might be visited.
namespace siconos::modeling {

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

namespace nonsmooth_laws {

struct Visitor {
  virtual ~Visitor() noexcept = default;
  // Note FP: use variant for nslaws ?
  // Pros: less virtual functions (no more class hierarchy in nslaws), can use std::visit
  // Cons: must use std::visit. Won't change anything at user level, but it will be a bit
  // more difficult to call functions at low level.

  virtual void visit(const siconos::modeling::ComplementarityConditionNSL &nslaw) {
    THROW_EXCEPTION("you must define a visit function in a derived class of NSLVisitor");
  };
  virtual void visit(const siconos::modeling::EqualityConditionNSL &nslaw) {
    THROW_EXCEPTION("you must define a visit function in a derived class of NSLVisitor");
  };
  virtual void visit(const siconos::modeling::MixedComplementarityConditionNSL &nslaw) {
    THROW_EXCEPTION("you must define a visit function in a derived class of NSLVisitor");
  };
  virtual void visit(const siconos::modeling::MultipleImpactNSL &nslaw) {
    THROW_EXCEPTION("you must define a visit function in a derived class of NSLVisitor");
  };
  virtual void visit(const siconos::modeling::NewtonImpactFrictionNSL &nslaw) {
    THROW_EXCEPTION("you must define a visit function in a derived class of NSLVisitor");
  };
  virtual void visit(const siconos::modeling::NewtonImpactNSL &nslaw) {
    THROW_EXCEPTION("you must define a visit function in a derived class of NSLVisitor");
  };
  virtual void visit(const siconos::modeling::NewtonImpactRollingFrictionNSL &nslaw) {
    THROW_EXCEPTION("you must define a visit function in a derived class of NSLVisitor");
  };
  virtual void visit(const siconos::modeling::NormalConeNSL &nslaw) {
    THROW_EXCEPTION("you must define a visit function in a derived class of NSLVisitor");
  };
  virtual void visit(const siconos::modeling::RelayNSL &nslaw) {
    THROW_EXCEPTION("you must define a visit function in a derived class of NSLVisitor");
  };
  virtual void visit(const siconos::modeling::FremondImpactFrictionNSL &nslaw) {
    THROW_EXCEPTION("you must define a visit function in a derived class of NSLVisitor");
  };
};

}  // namespace nonsmooth_laws
}  // namespace siconos::modeling

#endif /* NSLVisitor_hpp */
