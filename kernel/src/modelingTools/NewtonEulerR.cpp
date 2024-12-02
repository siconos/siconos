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
#include "PluginTypes.hpp"            // FPtr2 ...
#include "SiconosMatrixOp.hpp"        // setblock
#include "SiconosMatrixVectorOp.hpp"  // mat-vec prod
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"
#include "SimpleMatrix.hpp"
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::modeling::NewtonEulerR::initialize(Interaction& inter) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerR::initialize(Interaction& inter)\n");

  unsigned int ySize = inter.dimension();
  unsigned int xSize = inter.getSizeOfDS();
  unsigned int qSize = 7 * (xSize / 6);

  if (!_jachq)
    _jachq = std::make_shared<siconos::algebra::SimpleMatrix>(ySize, qSize);
  else {
    if (_jachq->size(0) == 0) {
      // if the matrix dim are null
      _jachq->resize(ySize, qSize);
    } else {
      assert((_jachq->size(1) == qSize && _jachq->size(0) == ySize) ||
             (printf("NewtonEuler::initializeWorkVectorsAndMatrices _jachq->size(1) = %d "
                     ",_qsize = %d , _jachq->size(0) = %d ,_ysize =%d \n",
                     _jachq->size(1), qSize, _jachq->size(0), ySize) &&
              false) ||
             ("NewtonEuler::initializeWorkVectorsAndMatrices inconsistent sizes between "
              "_jachq matrix and the interaction." &&
              false));
    }
  }

  DEBUG_EXPR(_jachq->display());

  if (!_jachqT) _jachqT = std::make_shared<siconos::algebra::SimpleMatrix>(ySize, xSize);

  if (!_T) {
    _T = std::make_shared<siconos::algebra::SimpleMatrix>(7, 6);
    _T->zero();
    _T->setValue(0, 0, 1.0);
    _T->setValue(1, 1, 1.0);
    _T->setValue(2, 2, 1.0);
  }
  DEBUG_EXPR(_jachqT->display());
  auto& DSlink = inter.linkToDSVariables();
  if (!_contactForce) {
    _contactForce = std::make_shared<siconos::algebra::SiconosVector>(
        DSlink[siconos::modeling::NewtonEulerR::p1]->size());
    _contactForce->zero();
  }
  DEBUG_END("siconos::modeling::NewtonEulerR::initialize(Interaction& inter)\n");
}

void siconos::modeling::NewtonEulerR::checkSize(Interaction& inter) {
  assert((_jachq->size(1) == 7 * (inter.getSizeOfDS() / 6) &&
          _jachq->size(0) == inter.dimension()) ||
         (printf("NewtonEuler::initializeWorkVectorsAndMatrices _jachq->size(1) = %d ,_qsize "
                 "= %d , _jachq->size(0) = %d ,_ysize =%d \n",
                 _jachq->size(1), 7 * (inter.getSizeOfDS() / 6), _jachq->size(0),
                 inter.dimension()) &&
          false) ||
         ("NewtonEuler::initializeWorkVectorsAndMatrices inconsistent sizes between _jachq "
          "matrix and the interaction." &&
          false));
}

void siconos::modeling::NewtonEulerR::setJachq(
    std::shared_ptr<siconos::algebra::SimpleMatrix> newJachq) {
  _jachq = newJachq;
}

void siconos::modeling::NewtonEulerR::setJachqPtr(
    std::shared_ptr<siconos::algebra::SimpleMatrix> newPtr) {
  _jachq = newPtr;
}

void siconos::modeling::NewtonEulerR::computeh(double time,
                                               const siconos::algebra::BlockVector& q0,
                                               siconos::algebra::SiconosVector& y) {
  siconos::algebra::prod(*_jachq, q0, y, true);
  if (_e) y += *_e;
}

void siconos::modeling::NewtonEulerR::computeOutput(double time, Interaction& inter,
                                                    unsigned int derivativeNumber) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerR::computeOutput(...)\n");
  DEBUG_PRINTF("with time = %f and derivativeNumber = %i starts\n", time, derivativeNumber);

  auto& DSlink = inter.linkToDSVariables();
  auto& y = *inter.y(derivativeNumber);
  auto& q = *DSlink[siconos::modeling::NewtonEulerR::q0];

  if (derivativeNumber == 0) {
    computeh(time, q, y);
  } else {
    /* \warning  V.A. 15/04/2016
     * We decide finally not to update the Jacobian there. To be discussed
     */
    // computeJachq(time, inter, DSlink[siconos::modeling::NewtonEulerR::q0]);
    // computeJachqT(inter, DSlink[siconos::modeling::NewtonEulerR::q0]);

    if (derivativeNumber == 1) {
      assert(_jachqT);
      assert(DSlink[siconos::modeling::NewtonEulerR::velocity]);
      DEBUG_EXPR(_jachqT->display(););
      DEBUG_EXPR((*DSlink[siconos::modeling::NewtonEulerR::velocity]).display(););

      siconos::algebra::prod(*_jachqT, *DSlink[siconos::modeling::NewtonEulerR::velocity], y);

      DEBUG_EXPR(y.display(););
    } else if (derivativeNumber == 2) {
      THROW_EXCEPTION(
          "Warning: we attempt to call siconos::modeling::NewtonEulerR::computeOutput(double "
          "time, Interaction& inter, InteractionProperties& interProp, unsigned int "
          "derivativeNumber) for derivativeNumber=2");
    } else
      THROW_EXCEPTION(
          "siconos::modeling::NewtonEulerR::computeOutput(double time, Interaction& inter, "
          "InteractionProperties& interProp, unsigned int derivativeNumber) derivativeNumber "
          "out of range or not yet implemented.");
  }
  DEBUG_END("siconos::modeling::NewtonEulerR::computeOutput(...)\n");
}

/** to compute p
 *  \param double : current time
 *  \Param unsigned int: "derivative" order of lambda used to compute input
 */
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

  if (level == 1) /* \warning : we assume that ContactForce is given by lambda[level] */
  {
    siconos::algebra::prod(lambda, *_jachqT, *_contactForce, true);

    DEBUG_PRINT("siconos::modeling::NewtonEulerR::computeInput contact force :\n");
    DEBUG_EXPR(_contactForce->display(););

    /*data is a pointer of memory associated to a dynamical system*/
    /** false because it consists in doing a sum*/
    siconos::algebra::prod(lambda, *_jachqT,
                           *DSlink[siconos::modeling::NewtonEulerR::p0 + level], false);

    DEBUG_EXPR(_jachqT->display(););
    DEBUG_EXPR(DSlink[siconos::modeling::NewtonEulerR::p0 + level]->display(););
    // DEBUG_EXPR(std::shared_ptr<siconos::algebra::SiconosVector> buffer =
    // std::make_shared<siconos::algebra::SiconosVector>(DSlink[siconos::modeling::NewtonEulerR::p0
    // + level]->size()));
    //            siconos::algebra::prod(lambda, *_jachqT, *buffer, true);
    //            std::cout << "added part to p" << buffer <<  std::endl;
    //            buffer->display(););
  }

  else if (level == 2) /* \warning : we assume that ContactForce is given by lambda[level] */
  {
    siconos::algebra::prod(lambda, *_jachqT, *_contactForce, true);
    DEBUG_EXPR(_contactForce->display(););

    /*data is a pointer of memory associated to a dynamical system*/
    /** false because it consists in doing a sum*/
    assert(DSlink[siconos::modeling::NewtonEulerR::p0 + level]);
    siconos::algebra::prod(lambda, *_jachqT,
                           *DSlink[siconos::modeling::NewtonEulerR::p0 + level], false);

    DEBUG_EXPR(_jachqT->display(););
    DEBUG_EXPR(DSlink[siconos::modeling::NewtonEulerR::p0 + level]->display(););
    DEBUG_EXPR(auto buffer = std::make_shared<siconos::algebra::SiconosVector>(
                   DSlink[siconos::modeling::NewtonEulerR::p0 + level]->size());
               siconos::algebra::prod(lambda, *_jachqT, *buffer, true);
               std::cout << "added part to p   " << buffer << std::endl; buffer->display(););
  } else if (level == 0) {
    siconos::algebra::prod(lambda, *_jachq,
                           *DSlink[siconos::modeling::NewtonEulerR::p0 + level], false);
  } else
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEulerR::computeInput(double time, Interaction& inter, "
        "InteractionProperties& interProp, unsigned int level)  not yet implemented for level "
        "> 1");
  DEBUG_END("siconos::modeling::NewtonEulerR::computeInput(...)\n");
}
/*It computes _jachqT=_jachq*T. Uploaded in the case of an unilateral constraint
 * (NewtonEuler3DR and NewtonEuler1DR)*/

void siconos::modeling::NewtonEulerR::computeJachqT(
    Interaction& inter, std::shared_ptr<siconos::algebra::BlockVector> q0) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerR::computeJachqT(Interaction& inter, "
      "std::shared_ptr<siconos::algebra::BlockVector> q0) \n");
  DEBUG_PRINTF("with inter =  %p\n", &inter);
  DEBUG_EXPR(inter.display());

  unsigned int k = 0;
  auto ySize = inter.dimension();
  auto auxBloc = std::make_shared<siconos::algebra::SimpleMatrix>(ySize, 7);
  auto auxBloc2 = std::make_shared<siconos::algebra::SimpleMatrix>(ySize, 6);
  std::vector<std::size_t> dimIndex(2);
  std::vector<std::size_t> startIndex(4);

  for (unsigned int i = 0; i < q0->numberOfBlocks(); i++) {
    auto q = (q0->getAllVect())[i];
    startIndex[0] = 0;
    startIndex[1] = 7 * k / 6;
    startIndex[2] = 0;
    startIndex[3] = 0;
    dimIndex[0] = ySize;
    dimIndex[1] = 7;
    siconos::algebra::setBlock(*_jachq, auxBloc, dimIndex, startIndex);

    computeT(q, _T);

    DEBUG_EXPR(q->display(););
    DEBUG_EXPR(_T->display());

    siconos::algebra::prod(*auxBloc, *_T, *auxBloc2);

    startIndex[0] = 0;
    startIndex[1] = 0;
    startIndex[2] = 0;
    startIndex[3] = k;
    dimIndex[0] = ySize;
    dimIndex[1] = 6;

    siconos::algebra::setBlock(*auxBloc2, _jachqT, dimIndex, startIndex);
    DEBUG_EXPR(_jachqT->display());

    k += 6;
  }
  DEBUG_END(
      "siconos::modeling::NewtonEulerR::computeJachqT(Interaction& inter, "
      "std::shared_ptr<siconos::algebra::BlockVector> q0) \n");
}

void siconos::modeling::NewtonEulerR::computeJach(double time, Interaction& inter) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerR::computeJachq(double time, Interaction& inter, ...) "
      "\n");
  DEBUG_PRINTF("with time =  %f\n", time);
  DEBUG_PRINTF("with inter =  %p\n", &inter);

  auto& DSlink = inter.linkToDSVariables();

  computeJachq(time, inter, DSlink[siconos::modeling::NewtonEulerR::q0]);
  computeJachqT(inter, DSlink[siconos::modeling::NewtonEulerR::q0]);

  // computeJachqDot(time, inter); // This is not needed here
  // computeDotJachq(time, inter);
  computeJachlambda(time, inter);
  DEBUG_END(
      "siconos::modeling::NewtonEulerR::computeJachq(double time, Interaction& inter, ...) "
      "\n");
}

void siconos::modeling::NewtonEulerR::computeDotJachq(
    double time, const siconos::algebra::BlockVector& workQ,
    siconos::algebra::BlockVector& workZ, const siconos::algebra::BlockVector& workQdot) {
  if (_dotjachq && _plugindotjacqh->fPtr) {
    auto zp = workZ.prepareVectorForPlugin();
    auto qp = workQ.prepareVectorForPlugin();
    auto qdotp = workQdot.prepareVectorForPlugin();
    ((siconos::plugins::FPtr2)(_plugindotjacqh->fPtr))(
        workQ.size(), qp->getArray(), workQdot.size(), qdotp->getArray(), &(*_dotjachq)(0, 0),
        zp->size(), &(*zp)(0));
    workZ = *zp;
  }
}

void siconos::modeling::NewtonEulerR::computeSecondOrderTimeDerivativeTerms(
    double time, Interaction& inter, std::shared_ptr<DynamicalSystem> ds1,
    std::shared_ptr<DynamicalSystem> ds2) {
  DEBUG_PRINT(
      "siconos::modeling::NewtonEulerR::computeSecondOrderTimeDerivativeTerms starts\n");

  auto& DSlink = inter.linkToDSVariables();
  // Compute the time derivative of the Jacobian
  if (!_dotjachq)  // lazy initialization
  {
    auto sizeY = inter.dimension();
    auto xSize = inter.getSizeOfDS();
    auto qSize = 7 * (xSize / 6);

    _dotjachq = std::make_shared<siconos::algebra::SimpleMatrix>(sizeY, qSize);
  }
  // Compute the product of the time derivative of the Jacobian with dotq
  // we assume that dotq is up to date !
  DEBUG_EXPR(DSlink[siconos::modeling::NewtonEulerR::dotq]->display(););
  computeDotJachq(time, *DSlink[siconos::modeling::NewtonEulerR::q0],
                  *DSlink[siconos::modeling::NewtonEulerR::z],
                  *DSlink[siconos::modeling::NewtonEulerR::dotq]);

  _secondOrderTimeDerivativeTerms =
      std::make_shared<siconos::algebra::SiconosVector>(_dotjachq->size(0));

  DEBUG_EXPR(_dotjachq->display(););

  siconos::algebra::prod(*_dotjachq, *DSlink[siconos::modeling::NewtonEulerR::dotq],
                         *_secondOrderTimeDerivativeTerms, true);

  DEBUG_EXPR(_secondOrderTimeDerivativeTerms->display());

  // Compute the product of jachq and Tdot --> jachqTdot

  unsigned int k = 0;
  auto ySize = inter.dimension();
  auto xSize = inter.getSizeOfDS();
  auto auxBloc = std::make_shared<siconos::algebra::SimpleMatrix>(ySize, 7);
  auto auxBloc2 = std::make_shared<siconos::algebra::SimpleMatrix>(ySize, 6);
  std::vector<std::size_t> dimIndex(2);
  std::vector<std::size_t> startIndex(4);

  auto jachqTdot = std::make_shared<siconos::algebra::SimpleMatrix>(ySize, xSize);
  bool endl = false;
  for (auto ds = ds1; !endl; ds = ds2) {
    endl = (ds == ds2);

    startIndex[0] = 0;
    startIndex[1] = 7 * k / 6;
    startIndex[2] = 0;
    startIndex[3] = 0;
    dimIndex[0] = ySize;
    dimIndex[1] = 7;
    siconos::algebra::setBlock(*_jachq, auxBloc, dimIndex, startIndex);

    NewtonEulerDS& d = *std::static_pointer_cast<NewtonEulerDS>(ds);
    d.computeTdot();
    auto& Tdot = *d.Tdot();

    DEBUG_EXPR(d.display());
    DEBUG_EXPR((d.Tdot())->display());

    siconos::algebra::prod(*auxBloc, Tdot, *auxBloc2);

    startIndex[0] = 0;
    startIndex[1] = 0;
    startIndex[2] = 0;
    startIndex[3] = k;
    dimIndex[0] = ySize;
    dimIndex[1] = 6;

    siconos::algebra::setBlock(*auxBloc2, jachqTdot, dimIndex, startIndex);
    DEBUG_EXPR(jachqTdot->display());

    k += ds->dimension();
  }

  // compute the product of jachqTdot and v
  siconos::algebra::SiconosVector workVelocity(
      *DSlink[siconos::modeling::NewtonEulerR::velocity]);
  DEBUG_EXPR(workVelocity.display(););
  siconos::algebra::prod(*jachqTdot, workVelocity, *_secondOrderTimeDerivativeTerms, false);
  DEBUG_EXPR(_secondOrderTimeDerivativeTerms->display());
  DEBUG_PRINT("siconos::modeling::NewtonEulerR::computeSecondOrderTimeDerivativeTerms ends\n");
}

void siconos::modeling::NewtonEulerR::accept(
    std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const {
  tourist->visit(*this);
}
