/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "Lagrangian2d2DR.hpp"
#include "Interaction.hpp"
#include "BlockVector.hpp"
#include "SimpleMatrix.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void Lagrangian2d2DR::initialize(Interaction& inter)
{
  //proj_with_q  _jachqProj.reset(new SimpleMatrix(_jachq->size(0),_jachq->size(1)));

  if((inter.getSizeOfDS() !=3) and (inter.getSizeOfDS() !=6))
  {
    THROW_EXCEPTION("Lagrangian2d2DR::initialize(Interaction& inter). The size of ds must of size 3");
  }
  unsigned int qSize = 3 * (inter.getSizeOfDS() / 3);
  _jachq = std::make_shared<SimpleMatrix>(2, qSize);
}

double Lagrangian2d2DR::distance() const
{
  DEBUG_BEGIN("Lagrangian2d2DR::distance(...)\n")
  SiconosVector dpc(*_Pc2 - *_Pc1);
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););
  DEBUG_EXPR(dpc.display(););
  DEBUG_END("Lagrangian2d2DR::distance(...)\n")
  return dpc.norm2() * (inner_prod(*_Nc, dpc) >= 0 ? -1 : 1);

}

void Lagrangian2d2DR::computeh(const BlockVector& q, BlockVector& z, SiconosVector& y)
{
  DEBUG_BEGIN("Lagrangian2d2DR::computeh(...)\n");
  DEBUG_EXPR(q.display());

  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););
  DEBUG_EXPR(_Nc->display(););

  LagrangianScleronomousR::computeh(q, z, y);
  y.setValue(0, distance());

  DEBUG_EXPR(y.display(););
  DEBUG_EXPR(display(););
  DEBUG_END("Lagrangian2d2DR::computeh(...)\n");
  //getchar();
}

void Lagrangian2d2DR::computeJachq(const BlockVector& q, BlockVector& z)
{
  DEBUG_BEGIN("Lagrangian2d2DR::computeJachq(Interaction& inter, SP::BlockVector q0 \n");

  double Nx = _Nc->getValue(0);
  double Ny = _Nc->getValue(1);
  double Px = _Pc1->getValue(0);
  double Py = _Pc1->getValue(1);
  double G1x = q.getValue(0);
  double G1y = q.getValue(1);


  /* construct tangent vector */
  double Tx = -Ny;
  double Ty =  Nx;


  double lever_arm_x = Px-G1x;
  double lever_arm_y = Py-G1y;
  DEBUG_PRINTF("N_x = %4.2e,\t N_ y = %4.2e\n", Nx, Ny);
  DEBUG_PRINTF("lever_arm_x = %4.2e,\t lever_arm_ y = %4.2e\n", lever_arm_x, lever_arm_y);

  // _jachq->setValue(0,0,Nx);
  // _jachq->setValue(0,1,Ny);
  // _jachq->setValue(0,2,lever_arm_x*Ny - lever_arm_y*Nx );

  // _jachq->setValue(1,0,Tx);
  // _jachq->setValue(1,1,Ty);
  // _jachq->setValue(1,2,lever_arm_x*Ty - lever_arm_y*Tx );

  double * array = &*_jachq->getArray();
  array[0] = Nx;
  array[2] = Ny;
  array[4] = lever_arm_x*Ny - lever_arm_y*Nx;

  array[1] = Tx;
  array[3] = Ty;
  array[5] = lever_arm_x*Ty - lever_arm_y*Tx;

  if(q.size() ==6)
  {
    DEBUG_PRINT("take into account second ds\n");
    double G2x = q.getValue(3);
    double G2y = q.getValue(4);
    lever_arm_x = Px-G2x;
    lever_arm_y = Py-G2y;

    DEBUG_PRINTF("lever_arm_x = %4.2e,\t lever_arm_ y = %4.2e\n", lever_arm_x, lever_arm_y);


    // _jachq->setValue(0,3,-Nx);
    // _jachq->setValue(0,4,-Ny);
    // _jachq->setValue(0,5,lever_arm_y * Nx - lever_arm_x*Ny);

    // _jachq->setValue(1,3,-Tx);
    // _jachq->setValue(1,4,-Ty);
    // _jachq->setValue(1,5,lever_arm_y * Tx - lever_arm_x*Ty);
    array[6] = -Nx;
    array[8] = -Ny;
    array[10] = lever_arm_y * Nx - lever_arm_x*Ny;

    array[7] = -Tx;
    array[9] = -Ty;
    array[11]= lever_arm_y*Tx - lever_arm_x*Ty;

  }
  DEBUG_EXPR(_jachq->display(););
  DEBUG_END("Lagrangian2d2DR::computeJachq(Interaction& inter, SP::BlockVector q0) \n");

}

void Lagrangian2d2DR::display() const
{
  LagrangianR::display();

  std::cout << " _Pc1 :" << std::endl;
  if(_Pc1)
    _Pc1->display();
  else
    std::cout << " nullptr :" << std::endl;

  std::cout << " _Pc2 :" << std::endl;
  if(_Pc2)
    _Pc2->display();
  else
    std::cout << " nullptr :" << std::endl;

  // std::cout << " _relPc1 :" << std::endl;
  // if(_relPc1)
  //   _relPc1->display();
  // else
  //   std::cout << " nullptr :" << std::endl;

  // std::cout << " _relPc2 :" << std::endl;
  // if(_relPc2)
  //   _relPc2->display();
  // else
  //   std::cout << " nullptr :" << std::endl;

  std::cout << " _Nc :" << std::endl;
  if(_Nc)
    _Nc->display();
  else
    std::cout << " nullptr :" << std::endl;
  // std::cout << " _relNc :" << std::endl;
  // if(_relNc)
  //   _relNc->display();
  // else
  //   std::cout << " nullptr :" << std::endl;

}



// void Lagrangian2d2DR::computeOutput(double time, Interaction& inter,  unsigned int derivativeNumber)
// {

//   DEBUG_PRINTF("Lagrangian2d2DR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber) with time = %f and derivativeNumber = %i\n", time, derivativeNumber);
//   VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
//   SiconosVector& y = *inter.y(derivativeNumber);
//   if(derivativeNumber == 0)
//   {
//     computeh(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z], y);
//   }
//   else
//   {
//     computeJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
//     if(derivativeNumber == 1)
//     {
//       assert(_jachq);

//       // direct prod to save time
//       //prod(*_jachq, *DSlink[LagrangianR::q1], y);

//       double * A = &*_jachq->getArray();
//       SP::BlockVector v = DSlink[LagrangianR::q1];
//       double *  v_ds_1 = v->vector(0)->getArray();

//       y(0)= A[0]* v_ds_1[0] + A[2]* v_ds_1[1]  + A[4]* v_ds_1[2];
//       y(1)= A[1]* v_ds_1[0] + A[3]* v_ds_1[1]  + A[5]* v_ds_1[2];

//       if (v ->numberOfBlocks() >1 )
//       {
//         double *  v_ds_2 = v->vector(1)->getArray();
//         y(0) += A[6]* v_ds_2[0] + A[8]* v_ds_2[1]  + A[10]* v_ds_2[2];
//         y(1) += A[7]* v_ds_2[0] + A[9]* v_ds_2[1]  + A[11]* v_ds_2[2];

//       }


//     // }
//     //   else
//     //   {
//     //     y(0)= A[0]* (*v)(0) + A[2]* (*v)(1) + A[4]* (*v)(2)
//     //         + A[6]* (*v)(3) + A[8]* (*v)(4) + A[10]* (*v)(5);
//     //     y(1)= A[1]* (*v)(0) + A[3]* (*v)(1) + A[5]* (*v)(2)
//     //         + A[7]* (*v)(3) + A[9]* (*v)(4) + A[11]* (*v)(5);

//     //   }
//       // for (unsigned int i =0; i < 2; i++)
//       // {
//       //   y(i)= A[i]* (*v)(0);
//       //   for (unsigned int j =1; j<v->size(); j++)
//       //   {
//       //     y(i) += A[i+j*2] * (*v)(j);
//       //   }
//       // }
//     }
//     else
//       THROW_EXCEPTION("Lagrangian2d2DR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber), index out of range");
//   }
// }





// void LagrangianScleronomousR::computeInput(double time, Interaction& inter, unsigned int level)
// {
//   DEBUG_BEGIN("void LagrangianScleronomousR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level) \n");

//   DEBUG_PRINTF("level = %i\n", level);
//   VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
//   computeJachq(*DSlink[LagrangianR::q0], *DSlink[LagrangianR::z]);
//   // get lambda of the concerned interaction
//   SiconosVector& lambda = *inter.lambda(level);
//   DEBUG_EXPR(lambda.display(););
//   DEBUG_EXPR(_jachq->display(););
//   // data[name] += trans(G) * lambda
//   //prod(lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);

//   double * A = &*_jachq->getArray();
//   SP::BlockVector v = DSlink[LagrangianR::q1];
//   int v_size= v->size();
//   for (unsigned int i =0; i < 2; i++)
//   {
//     y(i)= A[i]* (*v)(0);
//     for (unsigned int j =1; j<v->size(); j++)
//     {
//       y(i) += A[i+j*2] * (*v)(j);
//     }
//   }


//   DEBUG_EXPR(DSlink[LagrangianR::p0 + level]->display(););
//   DEBUG_END("void LagrangianScleronomousR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level) \n");
// }
