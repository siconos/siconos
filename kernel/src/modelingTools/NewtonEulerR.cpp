/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include <cstdio>

#include "NewtonEulerR.hpp"
#include "Interaction.hpp"
#include "NewtonEulerDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES

#include <debug.h>


#include <iostream>
// NewtonEulerR(): Relation(RELATION::NewtonEuler, RELATION::NonLinearR)
// {


// }

void NewtonEulerR::initialize(Interaction& inter)
{


  DEBUG_BEGIN("NewtonEulerR::initialize(Interaction& inter)\n");

  unsigned int ySize = inter.dimension();
  unsigned int xSize = inter.getSizeOfDS();
  unsigned int qSize = 7 * (xSize / 6);

  if (!_jachq)
    _jachq.reset(new SimpleMatrix(ySize, qSize));
  else
  {
    if (_jachq->size(0) == 0)
    {
      // if the matrix dim are null
      _jachq->resize(ySize, qSize);
    }
    else
    {
      assert((_jachq->size(1) == qSize && _jachq->size(0) == ySize) ||
             (printf("NewtonEuler::initializeWorkVectorsAndMatrices _jachq->size(1) = %d ,_qsize = %d , _jachq->size(0) = %d ,_ysize =%d \n", _jachq->size(1), qSize, _jachq->size(0), ySize) && false) ||
             ("NewtonEuler::initializeWorkVectorsAndMatrices inconsistent sizes between _jachq matrix and the interaction." && false));
    }
  }

  DEBUG_EXPR(_jachq->display());

  if (! _jachqT)
    _jachqT.reset(new SimpleMatrix(ySize, xSize));

  if (! _T)
  {
    _T.reset(new SimpleMatrix(7, 6));
    _T->zero();
    _T->setValue(0, 0, 1.0);
    _T->setValue(1, 1, 1.0);
    _T->setValue(2, 2, 1.0);
  }
  DEBUG_EXPR(_jachqT->display());
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  if (!_contactForce)
  {
    _contactForce.reset(new SiconosVector(DSlink[NewtonEulerR::p1]->size()));
    _contactForce->zero();
  }
  DEBUG_END("NewtonEulerR::initialize(Interaction& inter)\n");
}


void NewtonEulerR::checkSize(Interaction& inter)
{
  assert((_jachq->size(1) == 7 * (inter.getSizeOfDS() / 6) && _jachq->size(0) ==  inter.dimension()) ||
         (printf("NewtonEuler::initializeWorkVectorsAndMatrices _jachq->size(1) = %d ,_qsize = %d , _jachq->size(0) = %d ,_ysize =%d \n", _jachq->size(1), 7 * (inter.getSizeOfDS() / 6), _jachq->size(0),  inter.dimension()) && false) ||
         ("NewtonEuler::initializeWorkVectorsAndMatrices inconsistent sizes between _jachq matrix and the interaction." && false));
}



void NewtonEulerR::setJachq(SP::SimpleMatrix newJachq)
{
  _jachq = newJachq;
}

void NewtonEulerR::setJachqPtr(SP::SimpleMatrix newPtr)
{
  _jachq = newPtr ;
}



void NewtonEulerR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  prod(*_jachq, q0, y, true);
  if (_e)
    y += *_e;
}



void NewtonEulerR::computeOutput(double time, Interaction& inter, unsigned int derivativeNumber)
{

  DEBUG_BEGIN("NewtonEulerR::computeOutput(...)\n");
  DEBUG_PRINTF("with time = %f and derivativeNumber = %i starts\n", time, derivativeNumber);

  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  SiconosVector& y = *inter.y(derivativeNumber);
  BlockVector& q = *DSlink[NewtonEulerR::q0];


  if (derivativeNumber == 0)
  {
    computeh(time, q, y);
  }
  else
  {
    /* \warning  V.A. 15/04/2016
     * We decide finally not to update the Jacobian there. To be discussed
     */
    // computeJachq(time, inter, DSlink[NewtonEulerR::q0]);
    // computeJachqT(inter, DSlink[NewtonEulerR::q0]);

    if (derivativeNumber == 1)
    {
      assert(_jachqT);
      assert(DSlink[NewtonEulerR::velocity]);
      DEBUG_EXPR(_jachqT->display();); DEBUG_EXPR((*DSlink[NewtonEulerR::velocity]).display(););

      prod(*_jachqT, *DSlink[NewtonEulerR::velocity], y);

      DEBUG_EXPR(y.display(););
    }
    else if (derivativeNumber == 2)
    {
      RuntimeException::selfThrow("Warning: we attempt to call NewtonEulerR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber) for derivativeNumber=2");
    }
    else
      RuntimeException::selfThrow("NewtonEulerR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber) derivativeNumber out of range or not yet implemented.");
  }
  DEBUG_END("NewtonEulerR::computeOutput(...)\n");
}



/** to compute p
*  \param double : current time
*  \Param unsigned int: "derivative" order of lambda used to compute input
*/
void NewtonEulerR::computeInput(double time, Interaction& inter, unsigned int level)
{

  DEBUG_BEGIN("NewtonEulerR::computeInput(...)\n")
  DEBUG_PRINTF("with time = %f and level = %i starts\n", time, level);
  DEBUG_EXPR(printf("interaction %p\n",&inter););
  DEBUG_EXPR(inter.display(););
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();


  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);

  DEBUG_EXPR(lambda.display(););
  DEBUG_EXPR(DSlink[NewtonEulerR::p0 + level]->display(););

  if (level == 1) /* \warning : we assume that ContactForce is given by lambda[level] */
  {

    prod(lambda, *_jachqT, *_contactForce, true);

    DEBUG_PRINT("NewtonEulerR::computeInput contact force :\n");
    DEBUG_EXPR(_contactForce->display(););

    /*data is a pointer of memory associated to a dynamical system*/
    /** false because it consists in doing a sum*/
    prod(lambda, *_jachqT, *DSlink[NewtonEulerR::p0 + level], false);

    DEBUG_EXPR(_jachqT->display(););
    DEBUG_EXPR(DSlink[NewtonEulerR::p0 + level]->display(););
    // DEBUG_EXPR(SP::SiconosVector buffer(new SiconosVector(DSlink[NewtonEulerR::p0 + level]->size()));
    //            prod(lambda, *_jachqT, *buffer, true);
    //            std::cout << "added part to p" << buffer <<  std::endl;
    //            buffer->display(););

  }

  else if (level == 2) /* \warning : we assume that ContactForce is given by lambda[level] */
  {

    prod(lambda, *_jachqT, *_contactForce, true);
    DEBUG_EXPR(_contactForce->display(););

    /*data is a pointer of memory associated to a dynamical system*/
    /** false because it consists in doing a sum*/
    assert(DSlink[NewtonEulerR::p0 + level]);
    prod(lambda, *_jachqT, *DSlink[NewtonEulerR::p0 + level], false);

    DEBUG_EXPR(_jachqT->display(););
    DEBUG_EXPR(DSlink[NewtonEulerR::p0 + level]->display(););
    DEBUG_EXPR(SP::SiconosVector buffer(new SiconosVector(DSlink[NewtonEulerR::p0 + level]->size()));
               prod(lambda, *_jachqT, *buffer, true);
               std::cout << "added part to p   " << buffer <<  std::endl;
               buffer->display(););
  }
  else if (level == 0)
  {
    prod(lambda, *_jachq, *DSlink[NewtonEulerR::p0 + level], false);
  }
  else
    RuntimeException::selfThrow("NewtonEulerR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)  not yet implemented for level > 1");
  DEBUG_END("NewtonEulerR::computeInput(...)\n");
}
/*It computes _jachqT=_jachq*T. Uploaded in the case of an unilateral constraint (NewtonEuler3DR and NewtonEuler1DR)*/

void NewtonEulerR::computeJachqT(Interaction& inter, SP::BlockVector q0)
{
  DEBUG_BEGIN("NewtonEulerR::computeJachqT(Interaction& inter, SP::BlockVector q0) \n");
  DEBUG_PRINTF("with inter =  %p\n",&inter);
  DEBUG_EXPR(inter.display());

  unsigned int k = 0;
  unsigned int ySize = inter.dimension();
  SP::SimpleMatrix auxBloc(new SimpleMatrix(ySize, 7));
  SP::SimpleMatrix auxBloc2(new SimpleMatrix(ySize, 6));
  Index dimIndex(2);
  Index startIndex(4);

  for (unsigned int i =0 ; i < q0->numberOfBlocks()  ; i++)
  {
    SP::SiconosVector q = (q0->getAllVect())[i];
    startIndex[0] = 0;
    startIndex[1] = 7 * k / 6;
    startIndex[2] = 0;
    startIndex[3] = 0;
    dimIndex[0] = ySize;
    dimIndex[1] = 7;
    setBlock(_jachq, auxBloc, dimIndex, startIndex);

    computeT(q,_T);

    DEBUG_EXPR(q->display(););
    DEBUG_EXPR(_T->display());

    prod(*auxBloc, *_T, *auxBloc2);

    startIndex[0] = 0;
    startIndex[1] = 0;
    startIndex[2] = 0;
    startIndex[3] = k;
    dimIndex[0] = ySize;
    dimIndex[1] = 6;

    setBlock(auxBloc2, _jachqT, dimIndex, startIndex);
    DEBUG_EXPR(_jachqT->display());

    k += 6;
  }
  DEBUG_END("NewtonEulerR::computeJachqT(Interaction& inter, SP::BlockVector q0) \n");
}

void NewtonEulerR::computeJach(double time, Interaction& inter)
{
  DEBUG_BEGIN("NewtonEulerR::computeJachq(double time, Interaction& inter, ...) \n");
  DEBUG_PRINTF("with time =  %f\n",time);
  DEBUG_PRINTF("with inter =  %p\n",&inter);

  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();

  computeJachq(time, inter, DSlink[NewtonEulerR::q0]);
  computeJachqT(inter, DSlink[NewtonEulerR::q0]);

  //computeJachqDot(time, inter); // This is not needed here
  //computeDotJachq(time, inter);
  computeJachlambda(time, inter);
  DEBUG_END("NewtonEulerR::computeJachq(double time, Interaction& inter, ...) \n");
}

void NewtonEulerR::computeDotJachq(double time, BlockVector& workQ, BlockVector& workZ, BlockVector& workQdot)
{
  if (_dotjachq && _plugindotjacqh->fPtr)
    {
      ((FPtr2)(_plugindotjacqh->fPtr))(workQ.size(), &(workQ)(0), workQdot.size(), &(workQdot)(0), &(*_dotjachq)(0, 0), workZ.size(), &(workZ)(0));
    }
}

void  NewtonEulerR::computeSecondOrderTimeDerivativeTerms(double time, Interaction& inter, VectorOfBlockVectors& DSlink, SP::DynamicalSystem ds1, SP::DynamicalSystem ds2)
{
  DEBUG_PRINT("NewtonEulerR::computeSecondOrderTimeDerivativeTerms starts\n");

  // Compute the time derivative of the Jacobian
    if (!_dotjachq) // lazy initialization
  {
    unsigned int sizeY = inter.dimension();
    unsigned int xSize = inter.getSizeOfDS();
    unsigned int qSize = 7 * (xSize / 6);

    _dotjachq.reset(new SimpleMatrix(sizeY, qSize));
  }
  // Compute the product of the time derivative of the Jacobian with dotq
  BlockVector workQdot = *DSlink[NewtonEulerR::dotq]; // we assume that dotq is up to date !
  BlockVector workQ = *DSlink[NewtonEulerR::q0]; // we assume that dotq is up to date !
  BlockVector workZ = *DSlink[NewtonEulerR::z]; // we assume that dotq is up to date !
  DEBUG_EXPR(workQdot.display(););

  computeDotJachq(time, workQ, workZ, workQdot);

  _secondOrderTimeDerivativeTerms.reset(new SiconosVector(_dotjachq->size(0)));

  DEBUG_EXPR(_dotjachq->display(););

  prod(1.0,*_dotjachq, workQdot, *_secondOrderTimeDerivativeTerms, true);

  DEBUG_EXPR(_secondOrderTimeDerivativeTerms->display());

  // Compute the product of jachq and Tdot --> jachqTdot

  unsigned int k = 0;
  unsigned int ySize = inter.dimension();
  unsigned int xSize = inter.getSizeOfDS();
  SP::SimpleMatrix auxBloc(new SimpleMatrix(ySize, 7));
  SP::SimpleMatrix auxBloc2(new SimpleMatrix(ySize, 6));
  Index dimIndex(2);
  Index startIndex(4);

  SP::SimpleMatrix jachqTdot(new SimpleMatrix(ySize, xSize));
  bool endl = false;
  for (SP::DynamicalSystem ds = ds1; !endl; ds = ds2)
  {
    endl = (ds == ds2);

    startIndex[0] = 0;
    startIndex[1] = 7 * k / 6;
    startIndex[2] = 0;
    startIndex[3] = 0;
    dimIndex[0] = ySize;
    dimIndex[1] = 7;
    setBlock(_jachq, auxBloc, dimIndex, startIndex);

    NewtonEulerDS& d = *std11::static_pointer_cast<NewtonEulerDS> (ds);
    d.computeTdot();
    SimpleMatrix& Tdot = *d.Tdot();

    DEBUG_EXPR(d.display());
    DEBUG_EXPR((d.Tdot())->display());

    prod(*auxBloc, Tdot, *auxBloc2);

    startIndex[0] = 0;
    startIndex[1] = 0;
    startIndex[2] = 0;
    startIndex[3] = k;
    dimIndex[0] = ySize;
    dimIndex[1] = 6;

    setBlock(auxBloc2, jachqTdot, dimIndex, startIndex);
    DEBUG_EXPR(jachqTdot->display());

    k += ds->dimension();
  }

  // compute the product of jachqTdot and v
  SiconosVector workVelocity = *DSlink[NewtonEulerR::velocity];
  DEBUG_EXPR(workVelocity.display(););
  prod(1.0, *jachqTdot, workVelocity, *_secondOrderTimeDerivativeTerms, false);
  DEBUG_EXPR(_secondOrderTimeDerivativeTerms->display());
  DEBUG_PRINT("NewtonEulerR::computeSecondOrderTimeDerivativeTerms ends\n");

  *DSlink[NewtonEulerR::z] = workZ;
}
