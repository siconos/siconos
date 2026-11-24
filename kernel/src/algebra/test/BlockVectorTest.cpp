// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "BlockVectorTest.hpp"

#include "BlockVector.hpp"
#include "SiconosVector.hpp"

using namespace boost::numeric::ublas;
using BlockVector = siconos::algebra::BlockVector;

CPPUNIT_TEST_SUITE_REGISTRATION(BlockVectorTest);

void BlockVectorTest::setUp() {}

void BlockVectorTest::tearDown() {}

void BlockVectorTest::testDefaultConstructor() {
  std::cout << "==================================" << std::endl;
  std::cout << "=== BlockVector tests start ...=== " << std::endl;
  std::cout << "==================================" << std::endl;
  std::cout << "--> Test: default constructor" << std::endl;
  BlockVector x{};
  // Then insert some vectors
  auto w = std::make_shared<siconos::algebra::SiconosVector>(3);
  w->setConstant(2);
  auto z = std::make_shared<siconos::algebra::SiconosVector>(5);
  z->setConstant(3);
  x.insertPtr(w);
  x.insertPtr(z);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", x.size() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", x.numberOfBlocks() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", x.vector(0) == w, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", x.vector(1) == z, true);

  for (Eigen::Index i = 0; i < 3; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", x(i) == (*w)(i), true);
  for (Eigen::Index i = 3; i < x.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", x(i) == (*z)(i - 3), true);

  // setZero
  x.setZero();
  for (Eigen::Index i = 0; i < x.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", x(i) == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", (w->array() == 0).all(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", (z->array() == 0).all(), true);

  // setConstant
  x.setConstant(12.3);
  for (Eigen::Index i = 0; i < x.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", x(i) == 12.3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", (w->array() == 12.3).all(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", (z->array() == 12.3).all(), true);

  // Try print
  siconos::algebra::print(x);

  // fill

  double a[8] = {1., 2, 3, 4, 5, 6, 7, 8};
  x.fill(8, a);
  for (Eigen::Index i = 0; i < x.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", x(i) == i + 1, true);

  siconos::algebra::SiconosVector vec{8};
  vec.setConstant(4.);
  x.fill(vec);
  for (Eigen::Index i = 0; i < x.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", x(i) == 4., true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", (w->array() == 4.).all(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testDefaultConstructor : ", (z->array() == 4.).all(), true);

  std::cout << "✅--> Constructor test ended with success." << std::endl;
}

void BlockVectorTest::testConstructor2() {
  std::cout << "--> Test: constructor 2." << std::endl;
  BlockVector x{3, 4};
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", x.size() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", x.numberOfBlocks() == 3, true);
  for (const auto& v : x)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v->size() == 4, true);

  std::cout << "✅--> Test: constructor 2." << std::endl;
}

void BlockVectorTest::testConstructor3() {
  std::cout << "--> Test: constructor 3." << std::endl;
  BlockVector x{3};
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", x.size() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", x.numberOfBlocks() == 3, true);
  for (const auto& v : x)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v == nullptr, true);
  for (std::size_t nb = 0; nb < x.numberOfBlocks(); ++nb)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", x.vector(nb) == nullptr, true);

  auto w = std::make_shared<siconos::algebra::SiconosVector>(3);
  w->setConstant(3.6);

  x.setVectorPtr(1, w);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", x.vector(1) == w, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", x.size() == 3, true);
  auto z = std::make_shared<siconos::algebra::SiconosVector>(4);
  z->setConstant(2.3);
  x.setVectorPtr(0, z);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", x.vector(0) == z, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", x.vector(1) == w, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", x.size() == 7, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", x(0) == 2.3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", x(6) == 3.6, true);

  std::cout << "✅--> Constructor 3 test ended with success." << std::endl;
}

void BlockVectorTest::testNorm() {
  std::cout << "--> Test: norm." << std::endl;
  auto w = std::make_shared<siconos::algebra::SiconosVector>(3);
  auto z = std::make_shared<siconos::algebra::SiconosVector>(5);

  (*w)(0) = 1;
  (*w)(1) = 2;
  (*w)(2) = 3;

  (*z)(0) = 1;
  (*z)(1) = 2;
  (*z)(2) = 3;
  (*z)(3) = 4;
  (*z)(4) = 5;

  siconos::algebra::BlockVector xB{};
  xB.insertPtr(w);
  xB.insertPtr(z);

  double n2 = xB.norm();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testNorm : ", n2 - sqrt(69) < std::numeric_limits<double>::epsilon(), true);

  auto ni = siconos::algebra::normInf(xB);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNorm : ", ni == 5, true);
  std::cout << "✅--> norm test ended with success." << std::endl;
}

void BlockVectorTest::testAssignment() {
  std::cout << "--> Test: assignment." << std::endl;

  // Assignment
  siconos::algebra::BlockVector xB{};
  auto w0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  w0->setConstant(2);
  auto z0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  z0->setConstant(3);
  xB.insertPtr(w0);
  xB.insertPtr(z0);

  siconos::algebra::BlockVector yB{2, 3};
  yB = xB;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", yB.vector(0) != w0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", yB.vector(1) != z0, true);
  for (Eigen::Index i = 0; i < 3; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", yB(i) == (*w0)(i), true);
  for (Eigen::Index i = 3; i < yB.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", yB(i) == (*z0)(i - 3), true);

  yB *= 1.3;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (w0->array() == 2.).all(), true);

  // The blocks of the destination keep their memory: a block vector made of views on
  // other vectors updates them (used by the integrators for their work vectors).
  auto a0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  a0->setZero();
  siconos::algebra::BlockVector view{};
  view.insertPtr(a0);
  siconos::algebra::BlockVector src{1, 3};
  src.vector(0)->setConstant(5.);
  view = src;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment (view) : ", view.vector(0) == a0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment (view) : ", (a0->array() == 5.).all(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment (view) : ", view.vector(0) != src.vector(0),
                               true);
  std::cout << "✅--> assignment test ended with success." << std::endl;
}

void BlockVectorTest::testOperators() {
  std::cout << "--> Test: operators." << std::endl;

  siconos::algebra::BlockVector xB{};
  auto w0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  w0->setConstant(2.);
  auto z0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  z0->setConstant(2.);
  xB.insertPtr(w0);
  xB.insertPtr(z0);

  // +=
  siconos::algebra::BlockVector yB{2, 3};
  yB.setConstant(100.);

  yB += xB;
  for (Eigen::Index i = 0; i < yB.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators : ", yB(i) == 102, true);

  yB *= 2.;
  for (Eigen::Index i = 0; i < yB.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators : ", yB(i) == 204, true);

  siconos::algebra::SiconosVector vec{yB.size()};
  vec.setConstant(12.);
  yB -= vec;
  for (Eigen::Index i = 0; i < yB.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators : ", yB(i) == 192, true);
  std::cout << "✅--> operators test ended with success." << std::endl;
}

void BlockVectorTest::testToSiconosVector() {
  std::cout << "--> Test: testToSiconosVector." << std::endl;

  siconos::algebra::BlockVector xB{};
  auto w0 = std::make_shared<siconos::algebra::SiconosVector>(3);
  w0->setConstant(2.);
  auto z0 = std::make_shared<siconos::algebra::SiconosVector>(4);
  z0->setConstant(1.8);
  xB.insertPtr(w0);
  xB.insertPtr(z0);

  auto vec = xB.toSiconosVector();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testToSiconosVector : ", vec.size() == 7., true);

  for (Eigen::Index i = 0; i < 3; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testToSiconosVector : ", vec(i) == 2., true);
  for (Eigen::Index i = 3; i < vec.size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testToSiconosVector : ", vec(i) == 1.8, true);
  std::cout << "✅--> testToSiconosVector test ended with success." << std::endl;
}

void BlockVectorTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of BlockVector Tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
