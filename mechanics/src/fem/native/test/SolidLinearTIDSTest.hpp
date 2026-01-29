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
#ifndef __SolidLinearTIDSTest__
#define __SolidLinearTIDSTest__

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "Material.hpp"
#include "Mesh.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosSerialization.hpp"

class SolidLinearTIDSTest : public CppUnit::TestFixture
{

private:
  
  ACCEPT_SERIALIZATION(SolidLinearTIDSTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SolidLinearTIDSTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildSolidLinearTIDS1);
  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildSolidLinearTIDS1();

  // Members

  std::shared_ptr<siconos::mechanics::fem::Mesh> mesh;
  std::map<unsigned int, std::shared_ptr<siconos::mechanics::fem::Material>> materials;
  std::shared_ptr<siconos::algebra::SiconosMatrix> S,B;

public:
  void setUp();
  void tearDown();

};

#endif




