// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "MechanicsIOTest.hpp"

#include "Disk.hpp"
#include "MechanicsIO.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "SiconosMatrix.hpp"
// #include "SiconosKernel.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(MechanicsIOTest);

void MechanicsIOTest::setUp() {}

void MechanicsIOTest::tearDown() {}

void MechanicsIOTest::test1() {
  siconos::algebra::SiconosVector3 q;
  siconos::algebra::SiconosVector3 vel;
  q << 0., 1., 1.;
  vel << 0., 0., 10.;

  auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0, 10);

  auto ds1 = std::make_shared<siconos::collision::native::bodies::Disk>(1, 1, q, vel);
  auto ds2 = std::make_shared<siconos::collision::native::bodies::Disk>(2, 2, q, vel);

  nsds->insertDynamicalSystem(ds1);
  nsds->insertDynamicalSystem(ds2);

  siconos::io::MechanicsIO io;

  auto positions = io.positions(*nsds);
  auto velocities = io.velocities(*nsds);

  // ids
  CPPUNIT_ASSERT(positions(0, 0) == ds1->number());
  CPPUNIT_ASSERT(velocities(0, 0) == ds1->number());
  CPPUNIT_ASSERT(positions(0, 1) == ds2->number());
  CPPUNIT_ASSERT(velocities(0, 1) == ds2->number());

  CPPUNIT_ASSERT(positions.col(0).segment(1, 3) == q);
  CPPUNIT_ASSERT(velocities.col(0).segment(1, 3) == vel);
  CPPUNIT_ASSERT(positions.col(1).segment(1, 3) == q);
  CPPUNIT_ASSERT(velocities.col(1).segment(1, 3) == vel);
}
