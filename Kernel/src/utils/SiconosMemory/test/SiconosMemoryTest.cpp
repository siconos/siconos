c/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "SiconosMemoryTest.h"
#include "SimpleVector.h"
#include "BlockVector.h"
using namespace std;

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMemoryTest);

void SiconosMemoryTest::setUp()
{
  sizeMem = 3;
  unsigned int sizeVect = 3;

  vector<double> v(sizeVect);
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;
  vector<double> w(sizeVect);
  w[0] = 4;
  w[1] = 5;
  w[2] = 6;
  vector<double> z(sizeVect);
  z[0] = 7;
  z[1] = 8;
  z[2] = 9;
  q1.reset(new SimpleVector(v));
  q2.reset(new SimpleVector(w));
  q3.reset(new SimpleVector(z));
  c1.reset(new BlockVector());
  c2.reset(new BlockVector());

  (boost::static_pointer_cast<BlockVector>(c1))->insert(*q1);
  (boost::static_pointer_cast<BlockVector>(c1))->insert(*q2);
  (boost::static_pointer_cast<BlockVector>(c2))->insert(*q3);

  V1.push_back(q1);
  V1.push_back(q2);
  V2.push_back(c1);
  V2.push_back(c2);
  V3.push_back(q2);
  V3.push_back(q1);

}

void SiconosMemoryTest::tearDown()
{
}

// Constructor: data=memorySize
void SiconosMemoryTest::testBuildMemory1()
{

  SP::SiconosMemory tmp1(new SiconosMemory(4));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->getMemorySize() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp1->getVectorMemory().size() == 0, true);
  cout << "SiconosMemoryTest >>> testBuildMemory1 .............................. OK\n ";
}
// Constructor: copy of a std::vector of siconos vectors
void SiconosMemoryTest::testBuildMemory2()
{
  SP::SiconosMemory tmp1(new SiconosMemory(V1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : memorysize OK", tmp1->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : size vector OK", tmp1->getVectorMemory().size() == 2, true);

  SP::SiconosMemory tmp2(new SiconosMemory(V2));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : memorysize OK", tmp2->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : nbVectorsInMemory OK", tmp2->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory2 : size vector OK", tmp2->getVectorMemory().size() == 2, true);
  cout << "SiconosMemoryTest >>> testBuildMemory2 .............................. OK\n ";
}

// Constructor: copy of a std::vector of siconos vectors + memory size

void SiconosMemoryTest::testBuildMemory3()
{
  SP::SiconosMemory tmp1(new SiconosMemory(V1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : memorysize OK", tmp1->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : size vector OK", tmp1->getVectorMemory().size() == 2, true);

  SP::SiconosMemory tmp2(new SiconosMemory(V2));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : memorysize OK", tmp2->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : nbVectorsInMemory OK", tmp2->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory3 : size vector OK", tmp2->getVectorMemory().size() == 2, true);
  cout << "SiconosMemoryTest >>> testBuildMemory3 .............................. OK\n ";
}

// Copy constructor
void SiconosMemoryTest::testBuildMemory4()
{
  SP::SiconosMemory tmp1(new SiconosMemory(V1));

  SP::SiconosMemory tmp2(new SiconosMemory(*tmp1));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory4 : memorysize OK", tmp2->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory4 : nbVectorsInMemory OK", tmp2->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory4 : size vector OK", tmp2->getVectorMemory().size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory4 : size vector OK", *(tmp2->getVectorMemory())[0] == *(tmp1->getVectorMemory())[0] , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory4 : size vector OK", *(tmp2->getVectorMemory())[1] == *(tmp1->getVectorMemory())[1] , true);
  cout << "SiconosMemoryTest >>> testBuildMemory4 .............................. OK\n ";
}

// setVectorMemory
void SiconosMemoryTest::testSetVectorMemory()
{
  SP::SiconosMemory tmp1(new SiconosMemory(V1));

  tmp1->setVectorMemory(V2);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : memorysize OK", tmp1->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildMemory1 : size vector OK", tmp1->getVectorMemory().size() == 2, true);
  cout << "SiconosMemoryTest >>> testSetVectorMemory.............................. OK\n ";

}

// getSiconosVector
void SiconosMemoryTest::testGetSiconosVector()
{
  SP::SiconosMemory tmp1(new SiconosMemory(V1));
  SP::SiconosVector tmp = tmp1->getSiconosVector(0);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : *v1 size OK", tmp->size() == V1[0]->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVector : v1 values OK", *tmp == *(V1[0]), true);
  cout << "SiconosMemoryTest >>> testGetSiconosVector .......................... OK\n ";
}

// swap

void SiconosMemoryTest::testSwap()
{
  SP::SiconosMemory tmp1(new SiconosMemory(2));
  tmp1->swap((V1[0]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(tmp1->getSiconosVector(0)) == *(V1[0]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 1, true);
  tmp1->swap((V1[1]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(tmp1->getSiconosVector(0)) == *(V1[1]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(tmp1->getSiconosVector(1)) == *(V1[0]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  tmp1->swap((V1[0]));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(tmp1->getSiconosVector(0)) == *(V1[0]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : vector OK", *(tmp1->getSiconosVector(1)) == *(V1[1]), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == 2, true);
  cout << "SiconosMemoryTest >>> testSwap ...................................... OK\n ";
}


void SiconosMemoryTest::testOperatorEqual()
{
  SP::SiconosMemory tmp1(new SiconosMemory(V1));
  SP::SiconosMemory tmp2(new SiconosMemory(V3));

  *tmp2 = *tmp1;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : memorysize OK", tmp1->getMemorySize() ==  tmp2->getMemorySize(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : nbVectorsInMemory OK", tmp1->getNbVectorsInMemory() == tmp2->getNbVectorsInMemory(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : size vector OK", tmp1->getVectorMemory().size() == tmp2->getVectorMemory().size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : size vector OK", tmp1->getVectorMemory() == tmp2->getVectorMemory(), true);
  cout << "SiconosMemoryTest >>> testOperatorEqual ...................................... OK\n ";
}

void SiconosMemoryTest::testMemoryException()
{
  /** SiconosMemory::SiconosMemory(int memorySize, vector<SP::SiconosVector> V)
   * We test the case where memorySize < V.size
   */

  SiconosMemory M(1, V1);
}


void SiconosMemoryTest::testMemoryException1()
{
  /** SP::SiconosVector SiconosMemory::getSiconosVector(int index)
   * We test the case where index > memorySize
   */

  SP::SiconosMemory M(new SiconosMemory(V1));

  SP::SiconosVector v;
  v = M->getSiconosVector(V1.size() + 10);
}

