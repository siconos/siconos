// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __NewtonEulerRelationsTest__
#define __NewtonEulerRelationsTest__

#include <cppunit/extensions/HelperMacros.h>

/** Checks the Jacobians of the Newton-Euler relations, with one and two bodies.
 *
 *  Contact relations (NewtonEulerNormalVelocityR, NewtonEulerVelocityR,
 * NewtonEulerVelocityAngularVelocityR):
 *  - H (H): for twists (v_i, Omega_i) (v_i absolute, Omega_i in the body frame),
 *    H * [v1, Omega1, v2, Omega2] must be equal to the relative velocity at the contact
 *    point pc1, expressed in the contact frame [n, t1, t2]:
 *
 *      u = v1 + (R1 Omega1) x (pc1 - G1) - v2 - (R2 Omega2) x (pc1 - G2),
 *
 *    followed, for NewtonEulerVelocityAngularVelocityR, by the tangential components of the
 *    relative angular velocity R1 Omega1 - R2 Omega2.
 *  - J_h (G, first row): with pc1 - pc2 parallel to n, it must be equal to the gradient
 *    of the gap (finite differences), and J_h T v must be equal to the first row of H v.
 *
 *  NewtonEulerR (linear relation, constant J_h): the second-order terms
 *  Jdot qdot + J_h Tdot v, with the function for Jdot set after initialize.
 */
class NewtonEulerRelationsTest : public CppUnit::TestFixture {
 private:
  CPPUNIT_TEST_SUITE(NewtonEulerRelationsTest);

  CPPUNIT_TEST(testNewtonEuler1DRJacobian);
  CPPUNIT_TEST(testNewtonEuler3DRJacobian);
  CPPUNIT_TEST(testNewtonEulerVelocityAngularVelocityRJacobian);
  CPPUNIT_TEST(testNewtonEuler1DRGapGradient);
  CPPUNIT_TEST(testNewtonEuler3DRGapGradient);
  CPPUNIT_TEST(testNewtonEulerVelocityAngularVelocityRGapGradient);
  CPPUNIT_TEST(testNewtonEulerRSecondOrderTerms);
  CPPUNIT_TEST_SUITE_END();

  void testNewtonEuler1DRJacobian();
  void testNewtonEuler3DRJacobian();
  void testNewtonEulerVelocityAngularVelocityRJacobian();
  void testNewtonEuler1DRGapGradient();
  void testNewtonEuler3DRGapGradient();
  void testNewtonEulerVelocityAngularVelocityRGapGradient();
  void testNewtonEulerRSecondOrderTerms();

 public:
  void setUp();
  void tearDown();
};

#endif
