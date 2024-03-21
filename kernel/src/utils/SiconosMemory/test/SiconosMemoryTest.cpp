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
#include "SiconosMemoryTest.hpp"

#include <iostream>

#include "BlockVector.hpp"
#include "SiconosVector.hpp"

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMemoryTest);

using SiconosMemory = siconos::algebra::SiconosMemory;

constexpr auto sizeVect = 3;

void SiconosMemoryTest::setUp()
{
  _sizeMem = 3;

  std::vector<double> v(sizeVect);
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;
  std::vector<double> w(sizeVect);
  w[0] = 4;
  w[1] = 5;
  w[2] = 6;
  std::vector<double> z(sizeVect);
  z[0] = 7;
  z[1] = 8;
  z[2] = 9;
  q1 = std::make_shared<siconos::algebra::SiconosVector>(v.size());
  q2 = std::make_shared<siconos::algebra::SiconosVector>(w.size());
  q3 = std::make_shared<siconos::algebra::SiconosVector>(z.size());
  for (int i = 0; i < 3; i++)
  {
    (*q1)(i) = v[i];
    (*q2)(i) = w[i];
    (*q3)(i) = z[i];
  }
  
  c1 = std::make_shared<siconos::algebra::BlockVector>();
  c2 = std::make_shared<siconos::algebra::BlockVector>();

  c1->insertPtr(q1);
  c1->insertPtr(q2);
  c2->insertPtr(q3);

  V1 = std::make_shared<std::vector<siconos::algebra::SiconosVector>>();
  V2 = std::make_shared<std::vector<siconos::algebra::SiconosVector>>();
  V3 = std::make_shared<std::vector<siconos::algebra::SiconosVector>>();

  V1->push_back(*q1);
  V1->push_back(*q2);
  //  V2->push_back(*c1);
  //  V2->push_back(*c2);
  V3->push_back(*q2);
  V3->push_back(*q1);
}

void SiconosMemoryTest::tearDown() {}

// Constructor: data=memorySize
void SiconosMemoryTest::testBuildMemory1()
{
  std::cout << "=====================================" << std::endl;
  std::cout << "===  SiconosMemory tests start ...=== " << std::endl;
  std::cout << "=====================================" << std::endl;

  std::cout << "--> Test: constructor 0." << std::endl;
  auto tmp1 = std::make_shared<SiconosMemory>(4, sizeVect);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->size() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK",
                               tmp1->nbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp1->size() == 4, true);
  std::cout << "-->  testBuildMemory1 ended with success." << std::endl;
}

// setVectorMemory
void SiconosMemoryTest::testSetVectorMemory()
{
  std::cout << "--> Test: setVectorMemory." << std::endl;
  //  auto tmp1=std::make_shared<SiconosMemory(*V1);
  //
  //  tmp1->setVectorMemory(*V2);
  //
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->size() == 2,
  //  true); CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : _nbVectorsInMemory OK",
  //  tmp1->nbVectorsInMemory() == 2, true); CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 :
  //  size vector OK", tmp1->vectorMemory()->size() == 2, true);
  std::cout << "-->  setVectorMemory test ended with success." << std::endl;
}

// getSiconosVector
void SiconosMemoryTest::testGetSiconosVector()
{
  std::cout << "--> Test: getSiconosVector." << std::endl;
  auto tmp1 = std::make_shared<SiconosMemory>(2, sizeVect);
  tmp1->swap((*V1)[0]);
  tmp1->swap((*V1)[1]);
  const siconos::algebra::SiconosVector& tmp = tmp1->getSiconosVector(1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : *v1 size OK",
                               tmp.size() == (*V1)[0].size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : v1 values OK", tmp == ((*V1)[0]), true);
  std::cout << "-->  getSiconosVector test ended with success." << std::endl;
}

// swap

void SiconosMemoryTest::testSwap()
{
  std::cout << "--> Test: swap." << std::endl;
  auto tmp1 = std::make_shared<SiconosMemory>(2, sizeVect);
  tmp1->swap(((*V1)[0]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", tmp1->getSiconosVector(0) == ((*V1)[0]),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : _nbVectorsInMemory OK",
                               tmp1->nbVectorsInMemory() == 1, true);
  tmp1->swap(((*V1)[1]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", tmp1->getSiconosVector(0) == ((*V1)[1]),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", tmp1->getSiconosVector(1) == ((*V1)[0]),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : _nbVectorsInMemory OK",
                               tmp1->nbVectorsInMemory() == 2, true);
  tmp1->swap(((*V1)[0]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", tmp1->getSiconosVector(0) == ((*V1)[0]),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", tmp1->getSiconosVector(1) == ((*V1)[1]),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : _nbVectorsInMemory OK",
                               tmp1->nbVectorsInMemory() == 2, true);
  std::cout << "-->  swap test ended with success." << std::endl;
}

void SiconosMemoryTest::End()
{
  //   std::cout <<"======================================" <<std::endl;
  //   std::cout <<" ===== End of SiconosMemory Tests ===== " <<std::endl;
  //   std::cout <<"======================================" <<std::endl;
}
