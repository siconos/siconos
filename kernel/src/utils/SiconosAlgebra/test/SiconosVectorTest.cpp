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
#include "SiconosVectorTest.hpp"

#include "SiconosVectorIterator.hpp"

using namespace boost::numeric::ublas;
namespace ublas = boost::numeric::ublas;

using SiconosVector = siconos::algebra::SiconosVector;

CPPUNIT_TEST_SUITE_REGISTRATION(SiconosVectorTest);

void SiconosVectorTest::setUp() {
  tol = 1e-12;

  size = 5;
  size1 = 3;
  size2 = 2;

  // size = size1 + size2;

  ref = std::make_shared<SiconosVector>(size);
  for (unsigned int i = 0; i < size; ++i) (*ref)(i) = i;
  vq.resize(size, 1);
  for (unsigned int i = 0; i < size; i++) vq[i] = i + 1;

  // const vectors used for operators test (ex: x and y in z = x + y)
  // "B" in name for BlockVectors
  x = std::make_shared<SiconosVector>(vq);
  y = std::make_shared<SiconosVector>(vq);
  tmp1 = std::make_shared<SiconosVector>(size1);
  tmp2 = std::make_shared<SiconosVector>(size2);
  for (unsigned int i = 0; i < size1; ++i) (*tmp1)(i) = i;
  for (unsigned int i = 0; i < size2; ++i) (*tmp2)(i) = 100 * i;

  xB = std::make_shared<siconos::algebra::BlockVector>(tmp1, tmp2);
  yB = std::make_shared<siconos::algebra::BlockVector>(tmp2, tmp1);

  // vectors used as results
  z = std::make_shared<SiconosVector>(size);
  tmp3 = std::make_shared<SiconosVector>(size2);
  tmp4 = std::make_shared<SiconosVector>(size1);
  zB = std::make_shared<siconos::algebra::BlockVector>();
  zB->insertPtr(tmp3);
  zB->insertPtr(tmp4);
}

void SiconosVectorTest::tearDown() {}

void SiconosVectorTest::testSetBlock() {
  std::cout << "--> Test: setBlock." << std::endl;

  unsigned int sizeB = 9;
  auto subBlock = std::make_shared<SiconosVector>(sizeB);
  unsigned int pos = 1;
  // copy ref into subBlock
  siconos::algebra::setBlock(pos, *ref, *subBlock);

  for (unsigned int i = pos; i < pos + 5; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE(
        "test setBlock : ", fabs((*subBlock)(i) - (*ref)(i - pos)) < tol, true);
  std::cout << "--> setBlock test ended with success." << std::endl;
}

void SiconosVectorTest::testOrthoBaseFromVector() {
  /* test orthoBaseFromVector */

  siconos::algebra::SiconosVector3 n, t, s;
  siconos::algebra::SiconosVector3 base = {-1., 0., 1.};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        n << base(i), base(j), base(k);
        if ((i == 1 && j == 1) && k == 1) {
        } else {
          auto info = siconos::algebra::orthoBaseFromVector(n, t, s);
          std::cout << "n: " << n << "\n";
          std::cout << "t: " << t << "\n";
          std::cout << "s: " << s << "\n";
          std::cout << "n.s: " << n.dot(s) << "\n";
          std::cout << "n.s: " << n.dot(t) << "\n";
          std::cout << "t.s: " << t.dot(s) << "\n";
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector", info, true);
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector", n.dot(s) == 0., true);
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector", n.dot(t) == 0., true);
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector", t.dot(s) == 0., true);
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector",
                                       std::abs(n.norm() - 1.) < 1e-14, true);
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector",
                                       std::abs(t.norm() - 1.) < 1e-14, true);
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector",
                                       std::abs(s.norm() - 1.) < 1e-14, true);
        }
      }
    }
  }
}

void SiconosVectorTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of SiconosVector Tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
