// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __SolidLinearTIDSTest__
#define __SolidLinearTIDSTest__

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "Material.hpp"
#include "Mesh.hpp"
#include "SiconosMatrix.hpp"

class SolidLinearTIDSTest : public CppUnit::TestFixture {
 private:
  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SolidLinearTIDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildSolidLinearTIDS1);
  CPPUNIT_TEST(testPlaneStressCompliance);
  CPPUNIT_TEST(testLumpedMass);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildSolidLinearTIDS1();
  void testPlaneStressCompliance();
  void testLumpedMass();

  // Members

  std::shared_ptr<siconos::mechanics::fem::Mesh> mesh;
  siconos::algebra::SiconosSparseMatrix S{}, B{};
  siconos::mechanics::fem::Material mat{7800, 210e9, 1. / 3};

 public:
  void setUp();
  void tearDown();
};

#endif
