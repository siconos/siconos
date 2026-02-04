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

#include "RotationQuaternion.hpp"

using SiconosVector = siconos::algebra::SiconosVector;

CPPUNIT_TEST_SUITE_REGISTRATION(SiconosVectorTest);

void SiconosVectorTest::setUp() {}

void SiconosVectorTest::tearDown() {}

void SiconosVectorTest::testOrthoBaseFromVector() {
  /* test orthoBaseFromVector */

  siconos::algebra::SiconosVector3 n, t, s;
  siconos::algebra::SiconosVector3 base = {-1., 0., 1.};
  double tol = 1e-15;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        n << base(i), base(j), base(k);
        if ((i == 1 && j == 1) && k == 1) {
        } else {
          auto info = siconos::geometry::orthoBaseFromVector(n, t, s);
          // std::cout << "n: " << n << "\n";
          // std::cout << "t: " << t << "\n";
          // std::cout << "s: " << s << "\n";
          // std::cout << "n.s: " << n.dot(s) << "\n";
          // std::cout << "n.t: " << n.dot(t) << "\n";
          // std::cout << "t.s: " << t.dot(s) << "\n";
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector", info, true);
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector", n.dot(s) < tol, true);
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector", n.dot(t) < tol, true);
          CPPUNIT_ASSERT_EQUAL_MESSAGE("test orthoBaseFromVector", t.dot(s) < tol, true);
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
  std::cout << "✅ test orthoBaseFromVector passed.\n";
}

void SiconosVectorTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of SiconosVector Tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
