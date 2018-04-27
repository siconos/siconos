/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "SiconosVector.hpp"
#include "BlockVector.hpp"

using  std::cout;
using std::endl;

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMemoryTest);

unsigned int sizeVect = 3;

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
  q1.reset(new SiconosVector(v));
  q2.reset(new SiconosVector(w));
  q3.reset(new SiconosVector(z));
  c1.reset(new BlockVector());
  c2.reset(new BlockVector());

  c1->insert(*q1);
  c1->insert(*q2);
  c2->insert(*q3);

  V1.reset(new MemoryContainer);
  V2.reset(new MemoryContainer);
  V3.reset(new MemoryContainer);

  V1->push_back(*q1);
  V1->push_back(*q2);
  //  V2->push_back(*c1);
  //  V2->push_back(*c2);
  V3->push_back(*q2);
  V3->push_back(*q1);
}

void SiconosMemoryTest::tearDown()
{}

// Constructor: data=memorySize
void SiconosMemoryTest::testBuildMemory1()
{
  std::cout << "=====================================" <<std::endl;
  std::cout << "===  SiconosMemory tests start ...=== " <<std::endl;
  std::cout << "=====================================" <<std::endl;

  std::cout << "--> Test: constructor 0." <<std::endl;
  SP::SiconosMemory tmp1(new SiconosMemory(4, sizeVect));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->getMemorySize() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", tmp1->nbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp1->size() == 4, true);
  std::cout << "-->  testBuildMemory1 ended with success." <<std::endl;
}

// Constructor: copy of a std::vector of siconos vectors
void SiconosMemoryTest::testBuildMemory2()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::SiconosMemory tmp1(new SiconosMemory(*V1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : memorysize OK", tmp1->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : nbVectorsInMemory OK", tmp1->nbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : size vector OK", tmp1->size() == 2, true);

  //  SP::SiconosMemory tmp2(new SiconosMemory(*V2));
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : memorysize OK", tmp2->getMemorySize() == 2, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : nbVectorsInMemory OK", tmp2->nbVectorsInMemory() == 2, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : size vector OK", tmp2->vectorMemory()->size() == 2, true);
  std::cout << "-->  testBuildMemory2 ended with success." <<std::endl;
}

// Constructor: copy of a std::vector of siconos vectors + memory size

void SiconosMemoryTest::testBuildMemory3()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::SiconosMemory tmp1(new SiconosMemory(2, *V1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : memorysize OK", tmp1->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : nbVectorsInMemory OK", tmp1->nbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : size vector OK", tmp1->size() == 2, true);

  //  SP::SiconosMemory tmp2(new SiconosMemory(2,*V2));
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : memorysize OK", tmp2->getMemorySize() == 2, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : nbVectorsInMemory OK", tmp2->nbVectorsInMemory() == 2, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : size vector OK", tmp2->vectorMemory()->size() == 2, true);
  std::cout << "-->  testBuildMemory3 ended with success." <<std::endl;
}

// setVectorMemory
void SiconosMemoryTest::testSetVectorMemory()
{
  std::cout << "--> Test: setVectorMemory." <<std::endl;
  //  SP::SiconosMemory tmp1(new SiconosMemory(*V1));
  //
  //  tmp1->setVectorMemory(*V2);
  //
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->getMemorySize() == 2, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : _nbVectorsInMemory OK", tmp1->nbVectorsInMemory() == 2, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp1->vectorMemory()->size() == 2, true);
  std::cout << "-->  setVectorMemory test ended with success." <<std::endl;

}

// getSiconosVector
void SiconosMemoryTest::testGetSiconosVector()
{
  std::cout << "--> Test: getSiconosVector." <<std::endl;
  SP::SiconosMemory tmp1(new SiconosMemory(*V1));
  const SiconosVector& tmp = tmp1->getSiconosVector(0);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : *v1 size OK", tmp.size() == (*V1)[0].size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : v1 values OK", tmp == ((*V1)[0]), true);
  std::cout << "-->  getSiconosVector test ended with success." <<std::endl;
}

// swap

void SiconosMemoryTest::testSwap()
{
  std::cout << "--> Test: swap." <<std::endl;
  SP::SiconosMemory tmp1(new SiconosMemory(2, sizeVect));
  tmp1->swap(((*V1)[0]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", tmp1->getSiconosVector(0) == ((*V1)[0]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : _nbVectorsInMemory OK", tmp1->nbVectorsInMemory() == 1, true);
  tmp1->swap(((*V1)[1]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", tmp1->getSiconosVector(0) == ((*V1)[1]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", tmp1->getSiconosVector(1) == ((*V1)[0]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : _nbVectorsInMemory OK", tmp1->nbVectorsInMemory() == 2, true);
  tmp1->swap(((*V1)[0]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", tmp1->getSiconosVector(0) == ((*V1)[0]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", tmp1->getSiconosVector(1) == ((*V1)[1]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : _nbVectorsInMemory OK", tmp1->nbVectorsInMemory() == 2, true);
  std::cout << "-->  swap test ended with success." <<std::endl;
}

void SiconosMemoryTest::End()
{
  //   std::cout <<"======================================" <<std::endl;
  //   std::cout <<" ===== End of SiconosMemory Tests ===== " <<std::endl;
  //   std::cout <<"======================================" <<std::endl;
}


