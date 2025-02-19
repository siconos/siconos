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

#include "MBTB_JointR.hpp"

#include <boost/math/quaternion.hpp>

#include "NewtonEulerDS.hpp"
#include "NewtonEulerJointR.hpp"
#include "RotationQuaternion.hpp"  // For rewriteVectorFromBodyToAbsoluteFrame
#include "SiconosVectorOp.hpp"     // For cross_product
#include "SimpleMatrix.hpp"
#include "op3x3.h"  // For orthoBaseFromVector
// #define MBTB_JOINTR_DEBUG

siconos::mechanisms::MBTB_JointR::MBTB_JointR() {
  _M = std::make_shared<siconos::algebra::SiconosMatrix>(6, 6);
  _F = std::make_shared<siconos::algebra::SiconosVector>(6);
}

/*
 *
 *
 * ML_A_abs = M_FA_A + M_FB_A = F2/\BA
 * the system is :
 * FL=F1+F2
 * F1.e0 = F2.e0
 * ML_A_abs . ei = F2 /\ BA .ei for i=1,2
 *
 */

void siconos::mechanisms::MBTB_JointR::computeEquivalentForces() {
  if (!_G0C1) return;
  auto q = _ds1->q_read();
  double q1 = q(3);
  double q2 = q(4);
  double q3 = q(5);
  double q4 = q(6);
  ::boost::math::quaternion<double> quattrf(q1, q2, q3, q4);
  ::boost::math::quaternion<double> cquattrf(q1, -q2, -q3, -q4);
  ::boost::math::quaternion<double> quatbuff;
  siconos::algebra::SiconosVector3 FL = _jointR->contactForce().head(3);

  siconos::algebra::SiconosVector3 ML_G = _jointR->contactForce().tail(3);
  siconos::algebra::SiconosVector3 ML_G_abs;
  auto spML_G_abs = std::make_shared<siconos::algebra::SiconosVector3>();
  *spML_G_abs = ML_G;
  siconos::geometry::rewriteVectorFromBodyToAbsoluteFrame(*_ds1->q(), *spML_G_abs);
  ML_G_abs = *spML_G_abs;
#ifdef MBTB_JOINTR_DEBUG
  printf("siconos::mechanisms::MBTB_JointR::computeEquivalentForces Blambda\n");
  Blambda->display();
  printf(
      "siconos::mechanisms::MBTB_JointR::computeEquivalentForces ML_G_abs (in "
      "abs frame):\n");
  ML_G_abs.display();
#endif

  siconos::algebra::SiconosVector AB0(3);
  AB0 = (*_G0C2) - (*_G0C1);
  ::boost::math::quaternion<double> quatAB0(0, AB0(0), AB0(1), AB0(2));
  quatbuff = quattrf * quatAB0 * cquattrf;
  siconos::algebra::SiconosVector AB(3);
  AB(0) = quatbuff.R_component_2();
  AB(1) = quatbuff.R_component_3();
  AB(2) = quatbuff.R_component_4();
#ifdef MBTB_JOINTR_DEBUG
  printf("siconos::mechanisms::MBTB_JointR::computeEquivalentForces AB0:\n");
  AB0.display();
  printf("siconos::mechanisms::MBTB_JointR::computeEquivalentForces AB:\n");
  AB.display();
#endif

  ::boost::math::quaternion<double> quatG0C1(0, (*_G0C1)(0), (*_G0C1)(1), (*_G0C1)(2));
  quatbuff = quattrf * quatG0C1 * cquattrf;
  siconos::algebra::SiconosVector3 GA;
  GA(0) = quatbuff.R_component_2();
  GA(1) = quatbuff.R_component_3();
  GA(2) = quatbuff.R_component_4();
#ifdef MBTB_JOINTR_DEBUG
  printf("siconos::mechanisms::MBTB_JointR::computeEquivalentForces GA:\n");
  GA.display();
#endif
  siconos::algebra::SiconosVector ML_A_abs(3);
  siconos::algebra::SiconosVector3 GA_FL;
  GA_FL = GA.cross(FL);
  ML_A_abs = ML_G_abs - GA_FL;
#ifdef MBTB_JOINTR_DEBUG
  printf("siconos::mechanisms::MBTB_JointR::computeEquivalentForces GA_FL:\n");
  GA_FL.display();
  printf("siconos::mechanisms::MBTB_JointR::computeEquivalentForces ML_A_abs :\n");
  ML_A_abs.display();
#endif
  double normAB = 1.0 / AB.norm2();
  double e0_x = normAB * AB(0);
  double e0_y = normAB * AB(1);
  double e0_z = normAB * AB(2);
#ifdef MBTB_JOINTR_DEBUG
  printf(
      "siconos::mechanisms::MBTB_JointR::computeEquivalentForces e0_x . "
      "ML_A_abs (must be zero):%e\n",
      e0_x * ML_A_abs(0) + e0_y * ML_A_abs(1) + e0_z * ML_A_abs(2));
#endif
  double e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;
  int info =
      orthoBaseFromVector(&e0_x, &e0_y, &e0_z, &e1_x, &e1_y, &e1_z, &e2_x, &e2_y, &e2_z);
  if (info) {
    std::cout << "something wrong happened in  orthoBaseFromVector" << std::endl;
  }
  /*
   * ML_A_abs = M_FA_A + M_FB_A = F2/\BA
   * the system is :
   * FL=F1+F2
   * F1.e0 = F2.e0
   * ML_A_abs . ei = F2 /\ BA .ei for i=1,2
   *
   */
  /*Fill _M*/
  _M->setZero();
  // FL= F1+F2.
  (*_M)(0, 0) = 1;
  (*_M)(0, 3) = 1;
  (*_M)(1, 1) = 1;
  (*_M)(1, 4) = 1;
  (*_M)(2, 2) = 1;
  (*_M)(2, 5) = 1;
  /*F1.e0=F2.e0*/
  (*_M)(3, 0) = e0_x;
  (*_M)(3, 1) = e0_y;
  (*_M)(3, 2) = e0_z;
  (*_M)(3, 3) = -e0_x;
  (*_M)(3, 4) = -e0_y;
  (*_M)(3, 5) = -e0_z;

  /*
   *             F2x  BAx    F2y*BAz-F2z*BAy
   * FB/\BA . ei=F2y/\BAy.ei=F2z*BAx-F2x*BAz .ei = F2x*(BAy
   * *ei_z-BAz*ei_y)+F2y*(BAz*ei_x-BAx*ei_z)+F2z*(BAx*ei_y-BAy*ei_x) F2z  BAz
   * F2x*BAy-F2y*BAx
   *
   */
  (*_M)(4, 3) = -AB(1) * e1_z + AB(2) * e1_y;
  (*_M)(4, 4) = -AB(2) * e1_x + AB(0) * e1_z;
  (*_M)(4, 5) = -AB(0) * e1_y + AB(1) * e1_x;

  (*_M)(5, 3) = -AB(1) * e2_z + AB(2) * e2_y;
  (*_M)(5, 4) = -AB(2) * e2_x + AB(0) * e2_z;
  (*_M)(5, 5) = -AB(0) * e2_y + AB(1) * e2_x;
#ifdef MBTB_JOINTR_DEBUG
  printf(
      "siconos::mechanisms::MBTB_JointR::computeEquivalentForces the sytem M "
      "is:\n");
  _M->display();
#endif
  (*_F)(0) = FL(0);
  (*_F)(1) = FL(1);
  (*_F)(2) = FL(2);
  (*_F)(3) = 0;
  (*_F)(4) = e1_x * ML_A_abs(0) + e1_y * ML_A_abs(1) + e1_z * ML_A_abs(2);
  (*_F)(5) = e2_x * ML_A_abs(0) + e2_y * ML_A_abs(1) + e2_z * ML_A_abs(2);
#ifdef MBTB_JOINTR_DEBUG
  printf("siconos::mechanisms::MBTB_JointR::computeEquivalentForces Mx=b, b:");
  _F->display();
#endif

  /*Solve the system.*/
  try {
    auto luM = std::make_shared<Eigen::FullPivLU<siconos::algebra::SiconosMatrix>>(*_M);
    luM->solve(*_F);
#ifdef MBTB_JOINTR_DEBUG
    printf(
        "siconos::mechanisms::MBTB_JointR::computeEquivalentForces Forces "
        "equivalent:");
    _F->display();
    printf(
        "siconos::mechanisms::MBTB_JointR::computeEquivalentForces checking "
        "ML_G_abs = MF1_G+MF2_G:");
    siconos::algebra::SiconosVector Maux1(3), Maux2(3), F1(3), F2(3), GC1(3), GC2(3), dif(3);
    F1(0) = (*_F)(0);
    F1(1) = (*_F)(1);
    F1(2) = (*_F)(2);
    F2(0) = (*_F)(3);
    F2(1) = (*_F)(4);
    F2(2) = (*_F)(5);
    // FL=F1+F2
    dif = FL - F1 - F2;
    printf(
        "siconos::mechanisms::MBTB_JointR::computeEquivalentForces  "
        "FL-F1-F2(must be zero):\n");
    dif.display();
    // MF1_G=F1/\C1G
    Maux1 = GA.cross(F1);
    ::boost::math::quaternion<double> quatG0C2(0, (*_G0C2)(0), (*_G0C2)(1), (*_G0C2)(2));
    quatbuff = quattrf * quatG0C2 * cquattrf;
    siconos::algebra::SiconosVector GB(3);
    GB(0) = quatbuff.R_component_2();
    GB(1) = quatbuff.R_component_3();
    GB(2) = quatbuff.R_component_4();
    Maux2 = GA.cross(F2);
    dif = Maux1 + Maux2;
    printf(
        "siconos::mechanisms::MBTB_JointR::computeEquivalentForces  momentum "
        "(must be ML_G_abs):\n");
    dif.display();
    dif = Maux1 + Maux2 - ML_G_abs;
    printf(
        "siconos::mechanisms::MBTB_JointR::computeEquivalentForces  dif "
        "momentum(must be zero):\n");
    dif.display();
#endif
  } catch (const std::exception& e) {
    printf("MBTB_JointR: exception caught.\n");
  }
}
