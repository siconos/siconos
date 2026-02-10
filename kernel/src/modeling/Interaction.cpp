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
#include <memory>

#include "BlockVector.hpp"
#include "ComplementarityConditionNSL.hpp"
#include "DynamicalSystem.hpp"
#include "FirstOrderR.hpp"
#include "FremondImpactFrictionNSL.hpp"
#include "LagrangianR.hpp"
#include "NSLVisitor.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "NewtonImpactRollingFrictionNSL.hpp"
#include "Relation.hpp"
#include "RelayNSL.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMemory.hpp"
#include "SiconosVector.hpp"
//  #define DEBUG_BEGIN_END_ONLY
//   #define DEBUG_STDOUT
//   #define DEBUG_NOCOLOR
//   #define DEBUG_MESSAGES
#include "siconos_debug.h"

size_t siconos::modeling::Interaction::count_ = 0;

struct siconos::modeling::Interaction::SetLevels
    : public siconos::modeling::nonsmooth_laws::Visitor {
  /* we set the _lowerLevelForOutput, _upperLevelForOutput,
     _lowerLevelForOutput, _upperLevelForOutput
     w.r.t to the choice of the nslaw and the relation
  */
  using siconos::modeling::nonsmooth_laws::Visitor::visit;

  Interaction* interaction_{nullptr};

  SetLevels(Interaction* inter) : interaction_(inter) {};

  void visit(const ComplementarityConditionNSL& nslaw) override {
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

  void visit(const RelayNSL& nslaw) override {
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

  void visit(const NormalConeNSL& nslaw) override {
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

  void visit(const MixedComplementarityConditionNSL& nslaw) override {
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

  void visit(const EqualityConditionNSL& nslaw) override {
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
  void visit(const NewtonImpactNSL& nslaw) override {
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

  void visit(const NewtonImpactFrictionNSL& nslaw) override {
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
  void visit(const NewtonImpactRollingFrictionNSL& nslaw) override {
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

  void visit(const FremondImpactFrictionNSL& nslaw) override {
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

  void visit(const MultipleImpactNSL& nslaw) override {
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

  void visit(const MohrCoulombPlasticityNSL& nslaw) override {
    RelationType relationType = interaction_->relation()->getType();
    if (relationType == RelationType::Lagrangian ||
        relationType == RelationType::NewtonEuler) {
      interaction_->setLowerLevelForOutput(0);
      interaction_->setUpperLevelForOutput(0);

      interaction_->setLowerLevelForInput(0);
      interaction_->setUpperLevelForInput(0);
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
  auto nslawSize = _nslaw->size();

  for (auto i = _lowerLevelForOutput; i < _upperLevelForOutput + 1; i++) {
    _y[i] = std::make_shared<siconos::algebra::SiconosVector>(nslawSize);
    _y[i]->setZero();
  }

  for (auto i = _lowerLevelForInput; i < _upperLevelForInput + 1; i++) {
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

void siconos::modeling::Interaction::initialize_read_access_to_ds_variables(
    DynamicalSystem& ds1, DynamicalSystem& ds2) {
  std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars =
      read_access_to_ds_variables_;

  // The dynamical systems linked to the interaction (2 at most, ds2 may be equal to ds1).
  _relation->allocate_read_dynamical_systems_var_vectors(ds_vars, ds1, ds2);
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
  auto nslawSize = _nslaw->size();

  for (auto i = _lowerLevelForOutput; i < _upperLevelForOutput + 1; i++)
    _yMemory[i].setMemorySize(steps, nslawSize);

  for (auto i = _lowerLevelForInput; i < _upperLevelForInput + 1; i++) {
    DEBUG_PRINTF(
        "siconos::modeling::Interaction::initializeMemory(). "
        "_lambdaMemory[%i].setMemorySize()\n",
        i)
    _lambdaMemory[i].setMemorySize(steps, nslawSize);
  }

  DEBUG_END("siconos::modeling::Interaction::initializeMemory() \n");
}

void siconos::modeling::Interaction::resetAllLambda() {
  for (auto i = _lowerLevelForInput; i < _upperLevelForInput + 1; i++) {
    if (_lambda[i]) _lambda[i]->setZero();
  }
}

void siconos::modeling::Interaction::resetLambda(unsigned int level) {
  if (_lambda[level]) _lambda[level]->setZero();
}

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
  for (auto i = _lowerLevelForOutput; i < _upperLevelForOutput + 1; i++) {
    _yMemory[i].swap(*_y[i]);
  }
  for (auto i = _lowerLevelForInput; i < _upperLevelForInput + 1; i++) {
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

const siconos::algebra::ConstMapType siconos::modeling::Interaction::getLeftInteractionBlock()
    const {
  if (auto lagr = std::dynamic_pointer_cast<LagrangianR>(relation())) {
    return lagr->jacobianhOver_q();
  } else if (auto ner = std::dynamic_pointer_cast<NewtonEulerR>(relation())) {
    return ner->H_NE_prod_T();
  } else if (auto forel = std::dynamic_pointer_cast<FirstOrderR>(relation())) {
    return forel->jacobianhOver_state();
  } else {
    THROW_EXCEPTION(
        "siconos::modeling::Interaction::getLeftInteractionBlock, unknown relation type.");
  }
}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::modeling::Interaction::getLeftInteractionBlockForDS(
    siconos::algebra::Index pos, siconos::algebra::Index nslaw_size,
    siconos::algebra::Index ds_size) const {
  auto interactionBlock =
      std::make_shared<siconos::algebra::SiconosMatrix>(nslaw_size, ds_size);

  auto relationType = relation()->getType();
  if (relationType == RelationType::FirstOrder) {
    auto forel = std::dynamic_pointer_cast<FirstOrderR>(relation());
    if (forel->hasJacobianhOver_state()) {
      auto originalMatrix = forel->jacobianhOver_state();
      *interactionBlock = originalMatrix.block(0, pos, nslaw_size, ds_size);
    } else  // this should not happen? All first-order rel are supposed to have a
            // jacobianhOver_state
      interactionBlock->setZero();
  } else if (relationType == RelationType::Lagrangian) {
    auto lagr = std::dynamic_pointer_cast<LagrangianR>(relation());
    auto originalMatrix = lagr->jacobianhOver_q();
    *interactionBlock = originalMatrix.block(0, pos, nslaw_size, ds_size);

  } else if (relationType == RelationType::NewtonEuler) {
    auto newtonr = std::dynamic_pointer_cast<NewtonEulerR>(relation());
    auto originalMatrix = newtonr->H_NE_prod_T();
    *interactionBlock = originalMatrix.block(0, pos, nslaw_size, ds_size);
  } else
    THROW_EXCEPTION(
        "siconos::modeling::Interaction::getLeftInteractionBlockForDS, not yet implemented "
        "for relations of type " +
        std::to_string(
            static_cast<std::underlying_type<RelationSubType>::type>(relationType)));
  return interactionBlock;
}

void siconos::modeling::Interaction::getLeftInteractionBlockForDSProjectOnConstraints(
    siconos::algebra::Index pos,
    std::shared_ptr<siconos::algebra::SiconosMatrix> interactionBlock) const {
  DEBUG_PRINT(
      "siconos::modeling::Interaction::getLeftInteractionBlockForDSProjectOnConstraints()\n");
  DEBUG_PRINTF("pos = %i\n", pos);

  if (pos == 6) pos = pos + 1;
  auto neR = std::dynamic_pointer_cast<NewtonEulerR>(relation());
  assert(neR);
  // proj_with_q originalMatrix = r->jachqProj();
  auto originalMatrix = neR->H_NE();

  // copy sub-interactionBlock of originalMatrix into InteractionBlock
  *interactionBlock =
      originalMatrix.block(0, pos, interactionBlock->rows(), interactionBlock->cols());
}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::modeling::Interaction::getRightInteractionBlockForDS(
    siconos::algebra::Index pos, siconos::algebra::Index ds_size,
    siconos::algebra::Index nslaw_size) const {
  auto interactionBlock =
      std::make_shared<siconos::algebra::SiconosMatrix>(ds_size, nslaw_size);

  auto relationType = relation()->getType();

  if (relationType == RelationType::FirstOrder) {
    auto forel = std::dynamic_pointer_cast<FirstOrderR>(relation());
    assert(forel->hasJacobiangOver_lambda());
    auto originalMatrix = forel->jacobiangOver_lambda();
    *interactionBlock =
        originalMatrix.block(pos, 0, interactionBlock->rows(), interactionBlock->cols());
  } else  // (relationType == RelationType::Lagrangian ||
          //  relationType == RelationType::NewtonEuler) {
    THROW_EXCEPTION(
        "siconos::modeling::Interaction::getRightInteractionBlockForDS, unauthorized call "
        "for " +
        std::to_string(static_cast<std::underlying_type<RelationType>::type>(relationType)));
  return interactionBlock;
}

void siconos::modeling::Interaction::getExtraInteractionBlock(
    std::shared_ptr<siconos::algebra::SiconosMatrix> interactionBlock) const {
  // !!! Warning: we suppose that D is interactionBlock diagonal, ie that
  // there is no coupling between Interaction through D !!!  Any
  // coupling between relations through D must be taken into account
  // thanks to the nslaw (by "increasing" its dimension).

  if (relation()->hasJacobianhOver_lambda()) {
    auto originalMatrix = relation()->jacobianhOver_lambda();
    *interactionBlock = originalMatrix;  // copy!
  } else
    interactionBlock->setZero();
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
    std::cout << "| y[" << i << "] : \n";
    if (_y[i]) {
      if (_y[i]->size() >= 5) std::cout << "\n";
      siconos::algebra::print(*_y[i]);
    } else
      std::cout << "->nullptr\n";
  }
  for (unsigned int i = 0; i < _upperLevelForInput + 1; i++) {
    std::cout << "| lambda[" << i << "] :\n ";
    if (_lambda[i]) {
      if (_lambda[i]->size() >= 5) std::cout << "\n";
      siconos::algebra::print(*_lambda[i]);
    } else
      std::cout << "->nullptr\n";
  }
  if (!brief) {
    std::cout << "| _yMemory size: " << _yMemory.size() << "\n";

    for (unsigned int i = 0; i < _upperLevelForOutput + 1; i++) {
      std::cout << "| y_Memory[" << i << "] :\n ";
      siconos::algebra::print(_yMemory[i]);
    }
  }

  std::cout << "===================================\n";
}
