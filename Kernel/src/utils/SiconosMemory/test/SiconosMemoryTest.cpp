/* Siconos-Kernel, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "SiconosMemoryTest.hpp"
#include "SimpleVector.hpp"
#include "BlockVector.hpp"

using std::cout;
using std::endl;

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMemoryTest);

void SiconosMemoryTest::setUp()
{
  sizeMem = 3;
  unsigned int sizeVect = 3;

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
  q1.reset(new SimpleVector(v));
  q2.reset(new SimpleVector(w));
  q3.reset(new SimpleVector(z));
  c1.reset(new BlockVector());
  c2.reset(new BlockVector());

  c1->insert(*q1);
  c1->insert(*q2);
  c2->insert(*q3);

  V1.reset(new MemoryContainer);
  V2.reset(new MemoryContainer);
  V3.reset(new MemoryContainer);

  V1->push_back(q1);
  V1->push_back(q2);
  V2->push_back(c1);
  V2->push_back(c2);
  V3->push_back(q2);
  V3->push_back(q1);

}

void SiconosMemoryTest::tearDown()
{}

// Constructor: data=memorySize
void SiconosMemoryTest::testBuildMemory1()
{
  cout << "=====================================" << endl;
  cout << "===  SiconosMemory tests start ...=== " << endl;
  cout << "=====================================" << endl;

  cout << "--> Test: constructor 0." << endl;
  SP::SiconosMemory tmp1(new SiconosMemory(4));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->getMemorySize() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp1->vectorMemory()->size() == 0, true);
  cout << "-->  testBuildMemory1 ended with success." << endl;
}

// Constructor: copy of a std::vector of siconos vectors
void SiconosMemoryTest::testBuildMemory2()
{
  cout << "--> Test: constructor 1." << endl;
  SP::SiconosMemory tmp1(new SiconosMemory(*V1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : memorysize OK", tmp1->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : size vector OK", tmp1->vectorMemory()->size() == 2, true);

  SP::SiconosMemory tmp2(new SiconosMemory(*V2));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : memorysize OK", tmp2->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : nbVectorsInMemory OK", tmp2->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : size vector OK", tmp2->vectorMemory()->size() == 2, true);
  cout << "-->  testBuildMemory2 ended with success." << endl;
}

// Constructor: copy of a std::vector of siconos vectors + memory size

void SiconosMemoryTest::testBuildMemory3()
{
  cout << "--> Test: constructor 1." << endl;
  SP::SiconosMemory tmp1(new SiconosMemory(2, *V1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : memorysize OK", tmp1->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : size vector OK", tmp1->vectorMemory()->size() == 2, true);

  SP::SiconosMemory tmp2(new SiconosMemory(2, *V2));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : memorysize OK", tmp2->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : nbVectorsInMemory OK", tmp2->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : size vector OK", tmp2->vectorMemory()->size() == 2, true);
  cout << "-->  testBuildMemory3 ended with success." << endl;
}

// setVectorMemory
void SiconosMemoryTest::testSetVectorMemory()
{
  cout << "--> Test: setVectorMemory." << endl;
  SP::SiconosMemory tmp1(new SiconosMemory(*V1));

  tmp1->setVectorMemory(*V2);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp1->vectorMemory()->size() == 2, true);
  cout << "-->  setVectorMemory test ended with success." << endl;

}

// getSiconosVector
void SiconosMemoryTest::testGetSiconosVector()
{
  cout << "--> Test: getSiconosVector." << endl;
  SP::SiconosMemory tmp1(new SiconosMemory(*V1));
  SP::SiconosVector tmp = tmp1->getSiconosVector(0);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : *v1 size OK", tmp->size() == (*V1)[0]->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : v1 values OK", *tmp == *((*V1)[0]), true);
  cout << "-->  getSiconosVector test ended with success." << endl;
}

// swap

void SiconosMemoryTest::testSwap()
{
  cout << "--> Test: swap." << endl;
  SP::SiconosMemory tmp1(new SiconosMemory(2));
  tmp1->swap(((*V1)[0]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(tmp1->getSiconosVector(0)) == *((*V1)[0]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 1, true);
  tmp1->swap(((*V1)[1]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(tmp1->getSiconosVector(0)) == *((*V1)[1]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(tmp1->getSiconosVector(1)) == *((*V1)[0]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  tmp1->swap(((*V1)[0]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(tmp1->getSiconosVector(0)) == *((*V1)[0]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(tmp1->getSiconosVector(1)) == *((*V1)[1]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  cout << "-->  swap test ended with success." << endl;
}

void SiconosMemoryTest::End()
{
  //   cout <<"======================================" << endl;
  //   cout <<" ===== End of SiconosMemory Tests ===== " << endl;
  //   cout <<"======================================" << endl;
}


