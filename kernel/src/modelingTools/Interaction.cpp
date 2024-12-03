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

#include "Interaction.hpp"

#include <assert.h>

#include <iostream>

// #define DEBUG_BEGIN_END_ONLY
//  #define DEBUG_STDOUT
//  #define DEBUG_NOCOLOR
//  #define DEBUG_MESSAGES
#include "BlockVector.hpp"
#include "ComplementarityConditionNSL.hpp"
#include "DynamicalSystem.hpp"
#include "FirstOrderNonLinearDS.hpp"
#include "FirstOrderR.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianR.hpp"
// #include "FirstOrderLinearTIR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "NewtonImpactRollingFrictionNSL.hpp"
// #include "RelationType.hpp"
#include "RelayNSL.hpp"
#include "SiconosMatrix.hpp"
// include "SimulationGraphs.hpp"
#include "FremondImpactFrictionNSL.hpp"
#include "Relation.hpp"
#include "SiconosMatrixOp.hpp"  // For setBlock
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"
#include "siconos_debug.h"

// Test : the following line is allowed only from C++17.
#include <variant>

size_t siconos::modeling::Interaction::count_ = 0;

struct siconos::modeling::Interaction::SetLevels : public siconos::internal::SiconosVisitor {
  /* we set the _lowerLevelForOutput, _upperLevelForOutput,
     _lowerLevelForOutput, _upperLevelForOutput
     w.r.t to the choice of the nslaw and the relation
  */
  using SiconosVisitor::visit;

  Interaction* interaction_{nullptr};

  SetLevels(Interaction* inter) : interaction_(inter){};

  void visit(const ComplementarityConditionNSL& nslaw) const override {
    auto relationType = interaction_->relation()->getType();
    auto subType = interaction_->relation()->getSubType();

    if (relationType == RelationType::FirstOrder) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(0);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(0);
    } else if (relationType == RelationType::Lagrangian &&
               subType == RelationSubType::CompliantLinearTIR) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(0);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(0);
    } else {
      THROW_EXCEPTION(
          "siconos::modeling::Interaction::SetLevels::visit - unknown relation type for the "
          "nslaw ");
    };
  }

  void visit(const RelayNSL& nslaw) const override {
    RelationType relationType = interaction_->relation()->getType();
    if (relationType == RelationType::FirstOrder) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(0);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(0);
    } else if (relationType == RelationType::Lagrangian ||
               relationType == RelationType::NewtonEuler) {
      // For friction
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(1);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(1);
    } else {
      THROW_EXCEPTION(
          "siconos::modeling::Interaction::SetLevels::visit - unknown relation type for the "
          "nslaw ");
    };
  }

  void visit(const NormalConeNSL& nslaw) const override {
    RelationType relationType = interaction_->relation()->getType();
    if (relationType == RelationType::FirstOrder) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(0);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(0);
    } else {
      THROW_EXCEPTION(
          "siconos::modeling::Interaction::SetLevels::visit - unknown relation type for the "
          "nslaw ");
    };
  }

  void visit(const MixedComplementarityConditionNSL& nslaw) const override {
    RelationType relationType = interaction_->relation()->getType();
    if (relationType == RelationType::FirstOrder) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(0);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(0);
    } else {
      THROW_EXCEPTION(
          "siconos::modeling::Interaction::SetLevels::visit - unknown relation type for the "
          "nslaw ");
    };
  }

  void visit(const EqualityConditionNSL& nslaw) const override {
    RelationType relationType = interaction_->relation()->getType();
    if (relationType == RelationType::Lagrangian ||
        relationType == RelationType::NewtonEuler) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(1);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(1);
    } else if (relationType == RelationType::FirstOrder) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(0);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(0);
    } else {
      THROW_EXCEPTION(
          "siconos::modeling::Interaction::SetLevels::visit - unknown relation type for the "
          "nslaw ");
    };
  }
  void visit(const NewtonImpactNSL& nslaw) const override {
    RelationType relationType = interaction_->relation()->getType();
    if (relationType == RelationType::Lagrangian ||
        relationType == RelationType::NewtonEuler) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(1);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(1);
    } else {
      THROW_EXCEPTION(
          "siconos::modeling::Interaction::SetLevels::visit - unknown relation type for the "
          "nslaw ");
    }
  }

  void visit(const NewtonImpactFrictionNSL& nslaw) const override {
    RelationType relationType = interaction_->relation()->getType();
    if (relationType == RelationType::Lagrangian ||
        relationType == RelationType::NewtonEuler) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(1);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(1);
    } else {
      THROW_EXCEPTION(
          "siconos::modeling::Interaction::SetLevels::visit - unknown relation type for the "
          "nslaw ");
    }
  }
  void visit(const NewtonImpactRollingFrictionNSL& nslaw) const override {
    RelationType relationType = interaction_->relation()->getType();
    if (relationType == RelationType::Lagrangian ||
        relationType == RelationType::NewtonEuler) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(1);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(1);
    } else {
      THROW_EXCEPTION(
          "siconos::modeling::Interaction::SetLevels::visit - unknown relation type for the "
          "nslaw ");
    }
  }

  void visit(const FremondImpactFrictionNSL& nslaw) const override {
    auto relationType = interaction_->relation()->getType();
    if (relationType == RelationType::Lagrangian ||
        relationType == RelationType::NewtonEuler) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(1);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(1);

    } else {
      THROW_EXCEPTION("Interaction::_setLevels::visit - unknown relation type for the nslaw ");
    }
  }

  void visit(const MultipleImpactNSL& nslaw) const override {
    RelationType relationType = interaction_->relation()->getType();
    if (relationType == RelationType::Lagrangian ||
        relationType == RelationType::NewtonEuler) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(1);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(1);
    } else {
      THROW_EXCEPTION(
          "siconos::modeling::Interaction::SetLevels::visit - unknown relation type for the "
          "nslaw ");
    }
  }
};

void siconos::modeling::Interaction::reset() {
  // Check levels values and
  // resize all containers-like attributes according to these levels.

  // This function must be called at the first instanciation of
  // an interaction (in __init) and may be called by simulation and/or
  // OSI if levels are updated.

  assert(_upperLevelForOutput >= _lowerLevelForOutput);
  assert(_upperLevelForInput >= _lowerLevelForInput);

  // --  Memory allocation for y and lambda --
  // in order to simplify we size from 0 to _upperLevelForXXX
  _y.resize(_upperLevelForOutput + 1);
  _lambda.resize(_upperLevelForInput + 1);

  // get the dimension of the non smooth law, ie the size of an Interaction blocks (one per
  // relation)
  unsigned int nslawSize = _nslaw->size();

  for (unsigned int i = _lowerLevelForOutput; i < _upperLevelForOutput + 1; i++) {
    _y[i] = std::make_shared<siconos::algebra::SiconosVector>(nslawSize);
    _y[i]->setZero();
  }

  for (unsigned int i = _lowerLevelForInput; i < _upperLevelForInput + 1; i++) {
    _lambda[i] = std::make_shared<siconos::algebra::SiconosVector>(nslawSize);
  }
}

siconos::modeling::Interaction::Interaction(std::shared_ptr<NonSmoothLaw> NSL,
                                            std::shared_ptr<Relation> rel)
    : _number(count_++), _interactionSize(NSL->size()), _y(2), _nslaw(NSL), _relation(rel) {
  // -- Constructor --
  // i.e. what should be done when (and only there) the interaction
  // is instanciated.
  // Other operations (like levels review and y, lambda resizing)
  // occur in reset function, potentially called during
  // simulation phase (in OSI indeed).

  assert(_relation && "siconos::modeling::Interaction::__init failed, relation() == nullptr");
  assert(_nslaw && "siconos::modeling::Interaction::__inits, non smooth law == nullptr");

  // -- Set upper/lower levels, according to the nslaw --
  _nslaw->accept(*(std::make_shared<SetLevels>(this)));

  // Ensure consistency between interaction and nslaw sizes
  if (_interactionSize != _nslaw->size())
    THROW_EXCEPTION(
        "Interaction constructor - Nonsmooth law and relation are not consistent (sizes "
        "differ).");

  // Check levels and resize attributes (y, lambda ...) if needed.
  reset();
}

void siconos::modeling::Interaction::initializeLinkToDsVariables(DynamicalSystem& ds1,
                                                                 DynamicalSystem& ds2) {
  std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& DSlink = _linkToDSVariables;

  // The dynamical systems linked to the interaction (2 at most, ds2 may be equal to ds1).
  RelationType relationType = _relation->getType();

  if (relationType == RelationType::FirstOrder)
    __initDataFirstOrder(DSlink, ds1, ds2);

  else if (relationType == RelationType::Lagrangian)
    __initDataLagrangian(DSlink, ds1, ds2);

  else if (relationType == RelationType::NewtonEuler)
    __initDataNewtonEuler(DSlink, ds1, ds2);

  else
    THROW_EXCEPTION(
        "siconos::modeling::Interaction::initData unknown initialization procedure for \
        a relation of type: " +
        std::to_string(static_cast<std::underlying_type<RelationType>::type>(relationType)));

  _relation->initialize(*this);
}

// Initialize and InitializeMemory are separated in two functions
// since we need to know the relative degree to know
// "numberOfDerivatives", while numberOfRelations and the size of the
// non smooth law are required inputs to compute the relative degree.
void siconos::modeling::Interaction::initializeMemory(unsigned int steps) {
  DEBUG_BEGIN("siconos::modeling::Interaction::initializeMemory() \n");
  // Warning: this function is called from Simulation initialize,
  // since we need to know :
  // the levels _lowerLevelForOutput and _upperLevelForOutput to size Y
  // and the levels _lowerLevelForInput and _upperLevelForInput to size  Lambda.
  // this depends on many criteria (simulation type, osi type, ds type, nonsmooth type)
  // and they are computed in Simulation::computeLevelsForInputAndOutput

  // Memory allocation for y and lambda

  //  assert(_upperLevelForOutput >=0);
  assert(_upperLevelForOutput >= _lowerLevelForOutput);
  //  assert(_upperLevelForInput >=0);
  assert(_upperLevelForInput >= _lowerLevelForInput);

  _yMemory.resize(_upperLevelForOutput + 1);
  _lambdaMemory.resize(_upperLevelForInput + 1);
  unsigned int nslawSize = _nslaw->size();

  for (unsigned int i = _lowerLevelForOutput; i < _upperLevelForOutput + 1; i++)
    _yMemory[i].setMemorySize(steps, nslawSize);

  for (unsigned int i = _lowerLevelForInput; i < _upperLevelForInput + 1; i++) {
    DEBUG_PRINTF(
        "siconos::modeling::Interaction::initializeMemory(). "
        "_lambdaMemory[%i].setMemorySize()\n",
        i)
    _lambdaMemory[i].setMemorySize(steps, nslawSize);
  }

  DEBUG_END("siconos::modeling::Interaction::initializeMemory() \n");
}

void siconos::modeling::Interaction::resetAllLambda() {
  for (unsigned int i = _lowerLevelForInput; i < _upperLevelForInput + 1; i++) {
    if (_lambda[i]) _lambda[i]->zero();
  }
}

void siconos::modeling::Interaction::resetLambda(unsigned int level) {
  if (_lambda[level]) _lambda[level]->zero();
}

// It could be interesting to make Interaction a pure virtual class and to derive 3
// classes, one for each type of relation
void siconos::modeling::Interaction::__initDataFirstOrder(
    std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& DSlink, DynamicalSystem& ds1,
    DynamicalSystem& ds2) {
  DSlink.resize(FirstOrderR::DSlinkSize);
  DSlink[FirstOrderR::x] = std::make_shared<siconos::algebra::BlockVector>();
  DSlink[FirstOrderR::r] = std::make_shared<siconos::algebra::BlockVector>();
  DSlink[FirstOrderR::z] = std::make_shared<siconos::algebra::BlockVector>();
  auto relationSubType = _relation->getSubType();

  if (relationSubType != RelationSubType::LinearTIR) {
    // we need extra continuous memory vector
    // todo
  }

  __initDSDataFirstOrder(ds1, DSlink);
  if (&ds1 != &ds2) __initDSDataFirstOrder(ds2, DSlink);
}

void siconos::modeling::Interaction::__initDSDataFirstOrder(
    DynamicalSystem& ds, std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& DSlink) {
  // Put x/r ... of each DS into a block. (Pointers links, no copy!!)
  FirstOrderNonLinearDS& lds = static_cast<FirstOrderNonLinearDS&>(ds);
  DSlink[FirstOrderR::x]->insertPtr(lds.x());
  DSlink[FirstOrderR::r]->insertPtr(lds.r());
  DSlink[FirstOrderR::z]->insertPtr(lds.z());
}

void siconos::modeling::Interaction::__initDataLagrangian(
    std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& DSlink, DynamicalSystem& ds1,
    DynamicalSystem& ds2) {
  DEBUG_PRINT("siconos::modeling::Interaction::initDataLagrangian()\n");
  DSlink.resize(LagrangianR::DSlinkSize);

  // Default DSlink
  DSlink[LagrangianR::q0] = std::make_shared<siconos::algebra::BlockVector>();  // displacement
  DSlink[LagrangianR::q1] = std::make_shared<siconos::algebra::BlockVector>();  // velocity

  // auto relationSubType = _relation->getSubType();
  // if(relationSubType != LinearTIR)
  // {
  //   //we need extra continuous memory vector
  //   //todo
  // }

  __initDSDataLagrangian(ds1, DSlink);
  if (&ds1 != &ds2) __initDSDataLagrangian(ds2, DSlink);
}

void siconos::modeling::Interaction::__initDSDataLagrangian(
    DynamicalSystem& ds, std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& DSlink) {
  // check dynamical system type
  assert(dynamic_cast<LagrangianDS*>(&ds));

  auto& lds = static_cast<LagrangianDS&>(ds);

  // Put q, velocity of each DS into a block. (Pointers links, no copy!!)
  DSlink[LagrangianR::q0]->insertPtr(lds.q());
  DSlink[LagrangianR::q1]->insertPtr(lds.velocity());

  if (lds.acceleration()) {
    if (!DSlink[LagrangianR::q2])
      DSlink[LagrangianR::q2] =
          std::make_shared<siconos::algebra::BlockVector>();  // acceleration

    DSlink[LagrangianR::q2]->insertPtr(lds.acceleration());
  }

  if (lds.z()) {
    if (!DSlink[LagrangianR::z])
      DSlink[LagrangianR::z] = std::make_shared<siconos::algebra::BlockVector>();
    DSlink[LagrangianR::z]->insertPtr(lds.z());
  }
  for (unsigned int k = 0; k < 3; k++) {
    if (lds.p(k)) {
      if (!DSlink[LagrangianR::p0 + k])
        DSlink[LagrangianR::p0 + k] = std::make_shared<siconos::algebra::BlockVector>();
      DSlink[LagrangianR::p0 + k]->insertPtr(lds.p(k));
    }
  }
}

void siconos::modeling::Interaction::__initDataNewtonEuler(
    std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& DSlink, DynamicalSystem& ds1,
    DynamicalSystem& ds2) {
  DEBUG_BEGIN(
      "siconos::modeling::Interaction::initDataNewtonEuler(std::vector<std::shared_ptr<"
      "siconos::algebra::BlockVector>>& DSlink)\n");
  DSlink.resize(NewtonEulerR::DSlinkSize);
  // DSlink[NewtonEulerR::xfree] = std::make_shared<siconos::algebra::BlockVector>());
  DSlink[NewtonEulerR::q0] =
      std::make_shared<siconos::algebra::BlockVector>();  // displacement
  DSlink[NewtonEulerR::velocity] =
      std::make_shared<siconos::algebra::BlockVector>();                           // velocity
  DSlink[NewtonEulerR::dotq] = std::make_shared<siconos::algebra::BlockVector>();  // qdot
  //  data[NewtonEulerR::q2] = std::make_shared<siconos::algebra::BlockVector>(); //
  //  acceleration
  DSlink[NewtonEulerR::z] = std::make_shared<siconos::algebra::BlockVector>();  // z vector
  DSlink[NewtonEulerR::p0] = std::make_shared<siconos::algebra::BlockVector>();
  DSlink[NewtonEulerR::p1] = std::make_shared<siconos::algebra::BlockVector>();
  DSlink[NewtonEulerR::p2] = std::make_shared<siconos::algebra::BlockVector>();
  DEBUG_END(
      "siconos::modeling::Interaction::initDataNewtonEuler(std::vector<std::shared_ptr<"
      "siconos::algebra::BlockVector>>& DSlink)\n");
  __initDSDataNewtonEuler(ds1, DSlink);
  if (&ds1 != &ds2) __initDSDataNewtonEuler(ds2, DSlink);
}

void siconos::modeling::Interaction::__initDSDataNewtonEuler(
    DynamicalSystem& ds, std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& DSlink) {
  DEBUG_BEGIN(
      "siconos::modeling::Interaction::initDSDataNewtonEuler(DynamicalSystem& ds, "
      "std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& DSlink)\n");
  // check dynamical system type
  assert(
      (dynamic_cast<NewtonEulerDS*>(&ds)) &&
      "Interaction initDSData failed, not implemented for dynamical system of that type.\n");

  // convert vDS systems into NewtonEulerDS and put them in vLDS
  NewtonEulerDS& neds = static_cast<NewtonEulerDS&>(ds);
  // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
  DSlink[NewtonEulerR::q0]->insertPtr(neds.q());
  DSlink[NewtonEulerR::velocity]->insertPtr(neds.twist());
  //  DSlink[NewtonEulerR::deltaq]->insertPtr(neds.deltaq());
  DSlink[NewtonEulerR::dotq]->insertPtr(neds.dotq());
  //    data[NewtonEulerR::q2]->insertPtr( neds.acceleration());
  if (neds.p(0)) DSlink[NewtonEulerR::p0]->insertPtr(neds.p(0));
  if (neds.p(1)) DSlink[NewtonEulerR::p1]->insertPtr(neds.p(1));
  if (neds.p(2)) DSlink[NewtonEulerR::p2]->insertPtr(neds.p(2));

  DSlink[NewtonEulerR::z]->insertPtr(neds.z());
  DEBUG_END(
      "siconos::modeling::Interaction::initDSDataNewtonEuler(DynamicalSystem& ds, "
      "std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& DSlink)\n");
}
// --- GETTERS/SETTERS ---

void siconos::modeling::Interaction::setY(
    const std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>& newVector) {
  auto size = newVector.size();

  _y.clear();
  _y.resize(size);

  for (std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>::size_type i = 0;
       i < size; i++)
    _y[i] = std::make_shared<siconos::algebra::SiconosVector>(*(newVector[i]));  // -> copy !
}

void siconos::modeling::Interaction::setYPtr(
    const std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>& newVector) {
  _y.clear();

  // copy
  _y = newVector;  // smart ptr
}

void siconos::modeling::Interaction::setY(const unsigned int index,
                                          const siconos::algebra::SiconosVector& newY) {
  assert(_y.size() > index && "siconos::modeling::Interaction::setY, index out of range ");

  // set y[index]
  if (!_y[index]) {
    _y[index] = std::make_shared<siconos::algebra::SiconosVector>(newY);
  } else {
    assert(_y[index]->size() == newY.size() &&
           "siconos::modeling::Interaction::setY(index,newY), inconsistent sizes between "
           "y(index) and newY ");
    *(_y[index]) = newY;
  }
}

void siconos::modeling::Interaction::setYPtr(
    const unsigned int index, std::shared_ptr<siconos::algebra::SiconosVector> newY) {
  assert(_y.size() > index && "siconos::modeling::Interaction::setYPtr, index out of range");

  assert(newY->size() == _interactionSize &&
         "siconos::modeling::Interaction::setYPtr, interactionSize differs from newY vector "
         "size");

  _y[index] = newY;
}

void siconos::modeling::Interaction::setLambda(
    const std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>& newVector) {
  auto size = newVector.size();
  _lambda.clear();
  _lambda.resize(size);

  for (std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>::size_type i = 0;
       i < size; i++)
    _lambda[i] =
        std::make_shared<siconos::algebra::SiconosVector>(*(newVector[i]));  // -> copy !
}

void siconos::modeling::Interaction::setLambdaPtr(
    const std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>& newVector) {
  _lambda.clear();

  _lambda = newVector;  // smart ptr
}

void siconos::modeling::Interaction::setLambda(
    const unsigned int index, const siconos::algebra::SiconosVector& newLambda) {
  assert(_lambda.size() <= index &&
         "siconos::modeling::Interaction::setLambda, index out of range");

  // set lambda[index]
  if (!_lambda[index]) {
    _lambda[index] = std::make_shared<siconos::algebra::SiconosVector>(newLambda);
  } else {
    assert(_lambda[index]->size() == newLambda.size() &&
           "siconos::modeling::Interaction::setLambda(index,newLambda), inconsistent sizes "
           "between lambda(index) and newLambda");
    *(_lambda[index]) = newLambda;
  }
}

void siconos::modeling::Interaction::setLambdaPtr(
    const unsigned int index, std::shared_ptr<siconos::algebra::SiconosVector> newLambda) {
  assert(_lambda.size() > index &&
         "siconos::modeling::Interaction::setLambdaPtr, index out of range ");

  assert(newLambda->size() == _interactionSize &&
         "siconos::modeling::Interaction::setLambdaPtr, interactionSize differs from "
         "newLambda vector size ");

  _lambda[index] = newLambda;
}

const siconos::algebra::SiconosVector siconos::modeling::Interaction::getCopyOfy(
    const unsigned int i) const {
  assert(_y[i] && "_y[i]");
  return *(_y[i]);
}

siconos::algebra::SiconosMemory& siconos::modeling::Interaction::yMemory(unsigned int level) {
  return _yMemory[level];
}

const siconos::algebra::SiconosVector& siconos::modeling::Interaction::y_k(
    const unsigned int i) const {
  return _yMemory[i].getSiconosVector(0);
}

const siconos::algebra::SiconosVector siconos::modeling::Interaction::getLambda(
    const unsigned int i) const {
  assert(_lambda[i]);
  return *(_lambda[i]);
}

siconos::algebra::SiconosMemory& siconos::modeling::Interaction::lambdaMemory(
    unsigned int level) {
  return _lambdaMemory[level];
}

const siconos::algebra::SiconosVector& siconos::modeling::Interaction::lambda_k(
    const unsigned int i) const {
  return _lambdaMemory[i].getSiconosVector(0);
}

// --- OTHER FUNCTIONS ---

void siconos::modeling::Interaction::swapInMemory() {
  DEBUG_BEGIN("void siconos::modeling::Interaction::swapInMemory()\n");
  // i corresponds to the derivative number and j the relation number.
  for (unsigned int i = _lowerLevelForOutput; i < _upperLevelForOutput + 1; i++) {
    _yMemory[i].swap(*_y[i]);
  }
  for (unsigned int i = _lowerLevelForInput; i < _upperLevelForInput + 1; i++) {
    _lambdaMemory[i].swap(*_lambda[i]);
  }
  DEBUG_END("void siconos::modeling::Interaction::swapInMemory()\n");
}

void siconos::modeling::Interaction::computeOutput(double time,
                                                   unsigned int derivativeNumber) {
  DEBUG_BEGIN("siconos::modeling::Interaction::computeOutput(...)\n");
  DEBUG_PRINTF("time= %f\t", time);
  DEBUG_PRINTF("derivativeNumber= %i\n", derivativeNumber);
  relation()->computeOutput(time, *this, derivativeNumber);
  DEBUG_END("siconos::modeling::Interaction::computeOutput(...)\n");
}

void siconos::modeling::Interaction::computeInput(double time, unsigned int level) {
  DEBUG_BEGIN("siconos::modeling::Interaction::computeInput(...)\n");
  DEBUG_PRINTF("time= %f\t", time);
  DEBUG_PRINTF("level= %i\n", level);
  relation()->computeInput(time, *this, level);
  DEBUG_END("siconos::modeling::Interaction::computeInput(...)\n");
}

std::shared_ptr<siconos::algebra::MapType>
siconos::modeling::Interaction::getLeftInteractionBlock() const {
  auto relationType = relation()->getType();

  if (relationType == RelationType::Lagrangian) {
    std::shared_ptr<LagrangianR> r = std::static_pointer_cast<LagrangianR>(relation());
    return r->jachq();
  } else if (relationType == RelationType::NewtonEuler) {
    std::shared_ptr<NewtonEulerR> r = std::static_pointer_cast<NewtonEulerR>(relation());
    // return r->jachqT();
    return std::make_shared<siconos::algebra::MapType>( r->jachqT()->data(),  r->jachqT()->rows(),  r->jachqT()->cols()); // TODOSAM
  } else if (relationType == RelationType::FirstOrder) {
    std::shared_ptr<siconos::algebra::MapType> CMat =
        std::static_pointer_cast<FirstOrderR>(relation())->C(); // TODOSAM : not sure here
    auto relationSubType = relation()->getSubType();
    if (CMat)
      return CMat;
    else if (relationSubType != RelationSubType::LinearTIR)
      return std::make_shared<siconos::algebra::MapType>(_relationMatrices[FirstOrderR::mat_C]->data(), _relationMatrices[FirstOrderR::mat_C]->rows(), _relationMatrices[FirstOrderR::mat_C]->cols());
  }
  THROW_EXCEPTION(
      "siconos::modeling::Interaction::getLeftInteractionBlock, not yet implemented for "
      "relations of type " +
      std::to_string(static_cast<std::underlying_type<RelationType>::type>(relationType)));

  // return std::shared_ptr<siconos::algebra::SiconosMatrix>();
  return nullptr; //TODOSAM : handle this case
}

std::shared_ptr<siconos::algebra::MapType>
siconos::modeling::Interaction::getLeftInteractionBlockForDS(unsigned int pos, unsigned size,
                                                             unsigned int sizeDS) const {
  std::shared_ptr<siconos::algebra::MapType> originalMatrix;
  auto relationType = relation()->getType();
  if (relationType == RelationType::FirstOrder) {
    std::shared_ptr<siconos::algebra::MapType> CMat =
        std::static_pointer_cast<FirstOrderR>(relation())->C();
    auto relationSubType = relation()->getSubType();
    if (CMat)
      originalMatrix = CMat;
    else if (relationSubType != RelationSubType::LinearTIR)
      originalMatrix = std::make_shared<siconos::algebra::MapType>(_relationMatrices[FirstOrderR::mat_C]->data(), _relationMatrices[FirstOrderR::mat_C]->rows(), _relationMatrices[FirstOrderR::mat_C]->cols()); // TODOSAM
  } else if (relationType == RelationType::Lagrangian) {
    std::shared_ptr<LagrangianR> r = std::static_pointer_cast<LagrangianR>(relation());
    originalMatrix = r->jachq();
  } else if (relationType == RelationType::NewtonEuler) {
    std::shared_ptr<NewtonEulerR> r = std::static_pointer_cast<NewtonEulerR>(relation());
    originalMatrix = std::make_shared<siconos::algebra::MapType>( r->jachqT()->data(),  r->jachqT()->rows(),  r->jachqT()->cols()); // TODOSAM
  } else
    THROW_EXCEPTION(
        "siconos::modeling::Interaction::getLeftInteractionBlockForDS, not yet implemented "
        "for relations of type " +
        std::to_string(
            static_cast<std::underlying_type<RelationSubType>::type>(relationType)));

  auto InteractionBlock =
      std::make_shared<siconos::algebra::SiconosMatrix>(size, sizeDS);

  // copy sub-interactionBlock of originalMatrix into InteractionBlock
  // dim of the sub-interactionBlock
  std::vector<std::size_t> subDim(2);
  subDim[0] = InteractionBlock->size(0);
  subDim[1] = InteractionBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in InteractionBlock
  std::vector<std::size_t> subPos(4);
  subPos[0] = 0;  //_relativePosition;
  subPos[1] = pos;
  subPos[2] = 0;
  subPos[3] = 0;
  siconos::algebra::setBlock(*originalMatrix, InteractionBlock, subDim, subPos);
  return std::make_shared<siconos::algebra::MapType>(InteractionBlock->data(), InteractionBlock->rows(), InteractionBlock->cols()); // TODOSAM
}

void siconos::modeling::Interaction::getLeftInteractionBlockForDSProjectOnConstraints(
    unsigned int pos,
    std::shared_ptr<siconos::algebra::SiconosMatrix> InteractionBlock) const {
  DEBUG_PRINT(
      "siconos::modeling::Interaction::getLeftInteractionBlockForDSProjectOnConstraints("
      "unsigned int pos, std::shared_ptr<siconos::algebra::SiconosMatrix> InteractionBlock) "
      "\n");
  DEBUG_PRINTF("pos = %i\n", pos);

  if (pos == 6) pos = pos + 1;

  // Type::Siconos dsType = Type::value(*ds);
  // if (dsType != Type::NewtonEulerDS)
  //   THROW_EXCEPTION("siconos::modeling::Interaction::getLeftInteractionBlockForDSForProject-
  //   ds is not from NewtonEulerDS.");

  auto relationType = relation()->getType();
  if (relationType != RelationType::NewtonEuler)
    THROW_EXCEPTION(
        "siconos::modeling::Interaction::getLeftInteractionBlockForDSForProject- relation is "
        "not from NewtonEulerR.");

  std::shared_ptr<siconos::algebra::MapType> originalMatrix;
  std::shared_ptr<NewtonEulerR> r = std::static_pointer_cast<NewtonEulerR>(relation());
  // proj_with_q originalMatrix = r->jachqProj();
  originalMatrix = r->jachq();

  // copy sub-interactionBlock of originalMatrix into InteractionBlock
  // dim of the sub-interactionBlock
  std::vector<std::size_t> subDim(2);
  subDim[0] = InteractionBlock->size(0);
  subDim[1] = InteractionBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in InteractionBlock
  std::vector<std::size_t> subPos(4);
  subPos[0] = 0;  //_relativePosition;
  subPos[1] = pos;
  subPos[2] = 0;
  subPos[3] = 0;
  siconos::algebra::setBlock(*originalMatrix, InteractionBlock, subDim, subPos);
}

std::shared_ptr<siconos::algebra::MapType>
siconos::modeling::Interaction::getRightInteractionBlockForDS(unsigned int pos,
                                                              unsigned int sizeDS,
                                                              unsigned int size) const {
  std::shared_ptr<siconos::algebra::MapType>
      originalMatrix;  // Complete matrix, Relation member.
  auto relationType = relation()->getType();
  auto relationSubType = relation()->getSubType();

  if (relationType == RelationType::FirstOrder) {
    std::shared_ptr<siconos::algebra::MapType> BMat =
        std::static_pointer_cast<FirstOrderR>(relation())->B();
    if (BMat)
      originalMatrix = BMat;
    else if (relationSubType != RelationSubType::LinearTIR)
      originalMatrix = std::make_shared<siconos::algebra::MapType>(_relationMatrices[FirstOrderR::mat_B]->data(), _relationMatrices[FirstOrderR::mat_B]->rows(), _relationMatrices[FirstOrderR::mat_B]->cols()); // TODOSAM
    else
      THROW_EXCEPTION(
          "siconos::modeling::Interaction::getRightInteractionBlockForDS, FirstOrderLinearTIR "
          "relation but no B matrix found!");
  } else if (relationType == RelationType::Lagrangian ||
             relationType == RelationType::NewtonEuler) {
    THROW_EXCEPTION(
        "siconos::modeling::Interaction::getRightInteractionBlockForDS, call not permit " +
        std::to_string(static_cast<std::underlying_type<RelationType>::type>(relationType)));
  } else
    THROW_EXCEPTION(
        "siconos::modeling::Interaction::getRightInteractionBlockForDS, not yet implemented "
        "for relations of type " +
        std::to_string(
            static_cast<std::underlying_type<RelationSubType>::type>(relationType)));

  std::shared_ptr<siconos::algebra::SiconosMatrix> InteractionBlock =
      std::make_shared<siconos::algebra::SiconosMatrix>(sizeDS, size);

  if (!originalMatrix)
    THROW_EXCEPTION(
        "siconos::modeling::Interaction::getRightInteractionBlockForDS(DS, InteractionBlock, "
        "...): the right interactionBlock is a nullptr pointer (miss matrix B or H or "
        "gradients ...in relation ?)");

  // copy sub-interactionBlock of originalMatrix into InteractionBlock
  // dim of the sub-interactionBlock
  std::vector<std::size_t> subDim(2);
  subDim[0] = InteractionBlock->size(0);
  subDim[1] = InteractionBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in InteractionBlock
  std::vector<std::size_t> subPos(4);
  subPos[0] = pos;
  subPos[1] = 0;  //_relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  siconos::algebra::setBlock(*originalMatrix, InteractionBlock, subDim, subPos);
  return std::make_shared<siconos::algebra::MapType>(InteractionBlock->data(), InteractionBlock->rows(), InteractionBlock->cols()); // TODOSAM
}

void siconos::modeling::Interaction::getExtraInteractionBlock(
    std::shared_ptr<siconos::algebra::SiconosMatrix> InteractionBlock) const {
  // !!! Warning: we suppose that D is interactionBlock diagonal, ie that
  // there is no coupling between Interaction through D !!!  Any
  // coupling between relations through D must be taken into account
  // thanks to the nslaw (by "increasing" its dimension).

  RelationType relationType = relation()->getType();
  auto relationSubType = relation()->getSubType();
  std::shared_ptr<siconos::algebra::SiconosMatrix> D;

  if (relationType == RelationType::FirstOrder) {
    std::shared_ptr<siconos::algebra::SiconosMatrix> DMat =
        std::static_pointer_cast<FirstOrderR>(relation())->D();
    if (DMat)
      D = DMat;
    else if (relationSubType != RelationSubType::LinearTIR)
      D = _relationMatrices[FirstOrderR::mat_D];
  } else if (relationType == RelationType::Lagrangian) {
    // D = std::static_pointer_cast<LagrangianR>(relation())->jachlambda(); // TODOSAM : ?
  } else if (relationType == RelationType::NewtonEuler) {
    D = std::static_pointer_cast<NewtonEulerR>(relation())->jachlambda();
  } else
    THROW_EXCEPTION(
        "siconos::modeling::Interaction::getExtraInteractionBlockForDS, not yet implemented "
        "for relations of type " +
        std::to_string(
            static_cast<std::underlying_type<RelationSubType>::type>(relationType)));

  if (!D) {
    InteractionBlock->zero();
    return;  // ie no extra interactionBlock
  }

  *InteractionBlock = *D;
}
void siconos::modeling::Interaction::display(bool brief) const {
  std::cout << "======= Interaction display number " << _number << " =======\n";

  std::cout << "| lowerLevelForOutput : " << _lowerLevelForOutput << "\n";
  std::cout << "| upperLevelForOutput : " << _upperLevelForOutput << "\n";
  std::cout << "| lowerLevelForInput : " << _lowerLevelForInput << "\n";
  std::cout << "| upperLevelForInput : " << _upperLevelForInput << "\n";
  std::cout << "| interactionSize : " << _interactionSize << "\n";
  std::cout << "| _sizeOfDS : " << _sizeOfDS << "\n";

  std::cout << "| ";
  _relation->display();
  _nslaw->display();
  for (unsigned int i = 0; i < _upperLevelForOutput + 1; i++) {
    std::cout << "| y[" << i << "] : ";
    if (_y[i]) {
      if (_y[i]->size() >= 5) std::cout << "\n";
      _y[i]->display();
    } else
      std::cout << "->nullptr\n";
  }
  for (unsigned int i = 0; i < _upperLevelForInput + 1; i++) {
    std::cout << "| lambda[" << i << "] : ";
    if (_lambda[i]) {
      if (_lambda[i]->size() >= 5) std::cout << "\n";
      _lambda[i]->display();
    } else
      std::cout << "->nullptr\n";
  }
  if (!brief) {
    std::cout << "| _yMemory size: " << _yMemory.size() << "\n";
    ;
    for (unsigned int i = 0; i < _upperLevelForOutput + 1; i++) {
      std::cout << "| y_Memory[" << i << "] : ";
      _yMemory[i].display();
    }
  }

  std::cout << "===================================\n";
}
