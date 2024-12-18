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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "NewtonEulerR.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "NewtonEulerDS.hpp"  // computeT ...
#include "PluggedObject.hpp"
#include "PluginTypes.hpp"  // FPtr2 ...
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"  // setblock
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::modeling::NewtonEulerR::initialize(Interaction& inter) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerR::initialize(Interaction& inter)\n");

  unsigned int ySize = inter.dimension();
  unsigned int qSize = inter.getSizeOfDS();  // sum of considered DS sizes
  unsigned int Hcols = 7 * (qSize / 6);      // 7 * number of DS in the interaction

  // setHMatrix is the only way to set jacobianhOver_q. So if it's required, storage
  // is already allocated at this point.

  assert(jacobianhOver_q_view_);  // No relation without H ...

  // JacobianhOver_q_dot is optional. Allocation only if a user-defined function has been set.
  if (computejacobianhOver_q_dot_) {
    // True only if setComputeJacobianhOver_q_dotFunction has been called
    // Ensure that memory is properly allocated for H_dot
    if (!jacobianhOver_q_dot_) {
      jacobianhOver_q_dot_ = std::make_shared<siconos::algebra::SiconosMatrix>(ySize, Hcols);
    }
  }

  if (!jacobianhOver_q_prod_T_)
    jacobianhOver_q_prod_T_ = std::make_shared<siconos::algebra::SiconosMatrix>(ySize, qSize);

  // Allocate internal buffer, used to save T(q)
  if (!T_buffer_) {
    T_buffer_ = std::make_shared<siconos::algebra::SiconosMatrix>(7, 6);
    T_buffer_->setZero();
    T_buffer_->setValue(0, 0, 1.0);
    T_buffer_->setValue(1, 1, 1.0);
    T_buffer_->setValue(2, 2, 1.0);
  }

  if (!contactForce_) {
    auto& DSlink = inter.linkToDSVariables();
    contactForce_ = std::make_shared<siconos::algebra::SiconosVector>(
        DSlink[siconos::modeling::NewtonEulerR::p1]->size());
    contactForce_->setZero();
  }

  checkSize(inter);

  DEBUG_END("siconos::modeling::NewtonEulerR::initialize(Interaction& inter)\n");
}

void siconos::modeling::NewtonEulerR::checkSize(const Interaction &inter) const {
  unsigned int ySize = inter.dimension();
  unsigned int qSize = inter.getSizeOfDS();  // sum of considered DS sizes
  unsigned int Hcols = 7 * (qSize / 6);      // 7 * number of DS in the interaction

  if (jacobianhOver_q_view_) {
    assert(jacobianhOver_q_view_->rows() == ySize);
    assert(jacobianhOver_q_view_->cols() == Hcols);
  }

  if (jacobianhOver_q_dot_) {
    assert(jacobianhOver_q_dot_->rows() == ySize);
    assert(jacobianhOver_q_dot_->cols() == Hcols);
  }

  if (haseVector_) assert(eVector_view_->size() == ySize);

  if (jacobianhOver_q_prod_T_) {
    assert(jacobianhOver_q_prod_T_->rows() == ySize);
    assert(jacobianhOver_q_prod_T_->cols() == qSize);
  }
}

void siconos::modeling::NewtonEulerR::setConstantJacobianhOver_q(
    Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
  // Warning: mat. dims are not checked since we have no clue here regarding the size of the
  // Interaction.

  jacobianhOver_q_view_ = std::make_shared<siconos::algebra::MapType>(
      newValue.data(), newValue.rows(), newValue.cols());
}

void siconos::modeling::NewtonEulerR::setComputeJacobianhOver_q_dotFunction(
    const siconos::modeling::func_prototypes::FunctionBVBV_M& func) {
  // Ensure that memory is properly allocated for hmatrix_
  hasjacobianhOver_q_dot_ = true;
  computejacobianhOver_q_dot_ = func;
}

void siconos::modeling::NewtonEulerR::computeJacobianhOver_q_dot(
    const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot) {
  if (computejacobianhOver_q_dot_) {
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    assert(jacobianhOver_q_dot_);
    computejacobianhOver_q_dot_(q, qdot, *jacobianhOver_q_dot_);
  }
}

void siconos::modeling::NewtonEulerR::computeh(const siconos::algebra::BlockVector& q,
                                               Eigen::Ref<siconos::algebra::SiconosVector> y) {
  siconos::algebra::matrixBlockVector_prod(*jacobianhOver_q_view_, q, y, true);
  if (haseVector_) y += *eVector_view_;
}

void siconos::modeling::NewtonEulerR::computeOutput(double time, Interaction& inter,
                                                    unsigned int derivativeNumber) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerR::computeOutput(...)\n");
  DEBUG_PRINTF("with time = %f and derivativeNumber = %i starts\n", time, derivativeNumber);

  auto& DSlink = inter.linkToDSVariables();
  auto& y = *inter.y(derivativeNumber);
  auto& q = *DSlink[siconos::modeling::NewtonEulerR::q0];

  if (derivativeNumber == 0) {
    computeh(q, y);
  } else {
    /* \warning  V.A. 15/04/2016
     * We decide finally not to update the Jacobian there. To be discussed
     */
    // computeJacobianhOver_q(time, inter, DSlink[siconos::modeling::NewtonEulerR::q0]);
    // computeJacobianhOver_q_prod_T(inter, DSlink[siconos::modeling::NewtonEulerR::q0]);

    if (derivativeNumber == 1) {
      assert(jacobianhOver_q_prod_T_);
      assert(DSlink[siconos::modeling::NewtonEulerR::velocity]);
      DEBUG_EXPR(jacobianhOver_q_prod_T_->display(););
      DEBUG_EXPR((*DSlink[siconos::modeling::NewtonEulerR::velocity]).display(););

      siconos::algebra::matrixBlockVector_prod(
          *jacobianhOver_q_prod_T_, *DSlink[siconos::modeling::NewtonEulerR::velocity], y);

      DEBUG_EXPR(y.display(););
    } else if (derivativeNumber == 2) {
      THROW_EXCEPTION(
          "Call siconos::modeling::NewtonEulerR::computeOutput(..., derivativeNumber=2) is "
          "not allowed.");
    } else
      THROW_EXCEPTION(
          "siconos::modeling::NewtonEulerR::computeOutput(..., derivativeNumber) "
          "derivativeNumber out of range or not yet implemented.");
  }
  DEBUG_END("siconos::modeling::NewtonEulerR::computeOutput(...)\n");
}

void siconos::modeling::NewtonEulerR::computeInput(double time, Interaction& inter,
                                                   unsigned int level) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerR::computeInput(...)\n")
  DEBUG_PRINTF("with time = %f and level = %i starts\n", time, level);
  DEBUG_EXPR(printf("interaction %p\n", &inter););
  DEBUG_EXPR(inter.display(););
  auto& DSlink = inter.linkToDSVariables();

  // get lambda of the concerned interaction
  auto& lambda = *inter.lambda(level);

  DEBUG_EXPR(lambda.display(););
  DEBUG_EXPR(DSlink[siconos::modeling::NewtonEulerR::p0 + level]->display(););

  if (level == 1 ||
      level == 2) /* \warning : we assume that ContactForce is given by lambda[level] */
  {
    siconos::algebra::transposeMatrixVector_prod(lambda, *jacobianhOver_q_prod_T_,
                                                 *contactForce_, true);
    DEBUG_EXPR(contactForce_->display(););

    siconos::algebra::transposeMatrixVector_prod_toBlock(
        lambda, *jacobianhOver_q_prod_T_, *DSlink[siconos::modeling::NewtonEulerR::p0 + level],
        false);

    DEBUG_EXPR(jacobianhOver_q_prod_T_->display(););
    DEBUG_EXPR(DSlink[siconos::modeling::NewtonEulerR::p0 + level]->display(););
  } else if (level == 0) {
    siconos::algebra::transposeMatrixVector_prod_toBlock(
        lambda, *jacobianhOver_q_view_, *DSlink[siconos::modeling::NewtonEulerR::p0 + level],
        false);
  } else
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEulerR::computeInput(double time, Interaction& inter, "
        "InteractionProperties& interProp, unsigned int level)  not yet implemented for "
        "level "
        "> 1");
  DEBUG_END("siconos::modeling::NewtonEulerR::computeInput(...)\n");
}

void siconos::modeling::NewtonEulerR::computeJacobianhOver_q_prod_T(
    Interaction& inter, std::shared_ptr<siconos::algebra::BlockVector> q) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerR::computeJacobianhOver_q_prod_T(...) \n");
  DEBUG_PRINTF("with inter =  %p\n", &inter);
  DEBUG_EXPR(inter.display());

  unsigned int k = 0;
  auto ySize = inter.dimension();

  // For each qi from DS
  // compute H.T(qi)

  // H.T(q), first DS
  siconos::modeling::newton_euler::computeT(*q->vector(0), *T_buffer_);
  jacobianhOver_q_prod_T_->leftCols(6) = jacobianhOver_q_view_->leftCols(7) * *T_buffer_;
  // H.T(q), second DS if it exists
  if (inter.has2Bodies()) {
    siconos::modeling::newton_euler::computeT(*q->vector(1), *T_buffer_);
    jacobianhOver_q_prod_T_->rightCols(6) = jacobianhOver_q_view_->rightCols(7) * *T_buffer_;
  }
  DEBUG_END("siconos::modeling::NewtonEulerR::computeJacobianhOver_q_prod_T(...) \n");
}

void siconos::modeling::NewtonEulerR::computeJach(double time, Interaction& inter) {
  DEBUG_BEGIN("NewtonEulerR::computeJachq(double time, Interaction& inter, ...) \n");
  DEBUG_PRINTF("with time =  %f\n", time);
  DEBUG_PRINTF("with inter =  %p\n", &inter);

  auto& DSlink = inter.linkToDSVariables();

  computeJacobianhOver_q_prod_T(inter, DSlink[NewtonEulerR::q0]);

  // computeJachqDot(time, inter); // This is not needed here
  // computeDotJachq(time, inter);
  DEBUG_END("NewtonEulerR::computeJachq(double time, Interaction& inter, ...) \n");
}

void siconos::modeling::NewtonEulerR::computeSecondOrderTimeDerivativeTerms(
    double time, Interaction& inter, std::shared_ptr<DynamicalSystem> ds1,
    std::shared_ptr<DynamicalSystem> ds2) {
  DEBUG_PRINT(
      "siconos::modeling::NewtonEulerR::computeSecondOrderTimeDerivativeTerms starts\n");

  auto& DSlink = inter.linkToDSVariables();
  // Compute the time derivative of the Jacobian
  // and the product of the time derivative of the Jacobian with dotq
  // we assume that dotq is up to date !
  DEBUG_EXPR(DSlink[siconos::modeling::NewtonEulerR::dotq]->display();)
  assert(computejacobianhOver_q_dot_);
  computejacobianhOver_q_dot_(*DSlink[siconos::modeling::NewtonEulerR::q0],
                              *DSlink[siconos::modeling::NewtonEulerR::dotq],
                              *jacobianhOver_q_dot_);

  if (!secondOrderTimeDerivativeTerms_)
    secondOrderTimeDerivativeTerms_ =
        std::make_shared<siconos::algebra::SiconosVector>(jacobianhOver_q_dot_->size(0));

  DEBUG_EXPR(jacobianhOver_q_dot_->display(););

  siconos::algebra::matrixBlockVector_prod(*jacobianhOver_q_dot_,
                                           *DSlink[siconos::modeling::NewtonEulerR::dotq],
                                           *secondOrderTimeDerivativeTerms_, true);

  DEBUG_EXPR(_secondOrderTimeDerivativeTerms->display());

  // Compute the product of H and Tdot --> jachqTdot

  unsigned int k = 0;
  auto ySize = inter.dimension();
  auto qSize = inter.getSizeOfDS();
  auto auxBloc = std::make_shared<siconos::algebra::SiconosMatrix>(ySize, 7);

  std::vector<std::size_t> dimIndex(2);
  std::vector<std::size_t> startIndex(4);

  // Temp buffer to save H.Tdot
  siconos::algebra::SiconosMatrix jacobianhOver_q_prod_Tdot{ySize, qSize};
  auto neds = std::dynamic_pointer_cast<NewtonEulerDS>(ds1);
  neds->computeTdot();
  // H.Tdot(q)
  jacobianhOver_q_prod_Tdot.leftCols(6) = jacobianhOver_q_view_->leftCols(7) * neds->Tdot();
  if (ds2) {  // For second ds, if it exists
    neds = std::dynamic_pointer_cast<NewtonEulerDS>(ds2);
    neds->computeTdot();
    jacobianhOver_q_prod_Tdot.rightCols(6) =
        jacobianhOver_q_view_->rightCols(7) * neds->Tdot();
  }

  // compute the product of jachqTdot and v
  siconos::algebra::matrixBlockVector_prod(jacobianhOver_q_prod_Tdot,
                                           *DSlink[siconos::modeling::NewtonEulerR::velocity],
                                           *secondOrderTimeDerivativeTerms_, false);

  DEBUG_EXPR(_secondOrderTimeDerivativeTerms->display());
  DEBUG_PRINT("siconos::modeling::NewtonEulerR::computeSecondOrderTimeDerivativeTerms ends\n");
}

void siconos::modeling::NewtonEulerR::accept(
    std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const {
  tourist->visit(*this);
}
