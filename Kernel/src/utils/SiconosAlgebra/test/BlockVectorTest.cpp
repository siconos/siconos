/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

#include "BlockVectorTest.h"
#include "BlockVector.h"
#include "SimpleMatrix.h"
#include "SimpleVector.h"

using namespace std;
using namespace boost::numeric::ublas;

CPPUNIT_TEST_SUITE_REGISTRATION(BlockVectorTest);


void BlockVectorTest::setUp()
{
  tol = 1e-14;
  ref = new SimpleVector(5);
  for (unsigned int i = 0; i < 5; ++i)
    (*ref)(i) = i;

  vq.resize(5, 1);
  for (unsigned int i = 0; i < 5; i++)
    vq[i] = i + 1;

  dv = new DenseVect(3);
  (*dv)(0) = 1;
  (*dv)(1) = 2;
  (*dv)(2) = 3;
  sv = new SparseVect(5);
  (*sv)(1) = 22;
}

void BlockVectorTest::tearDown()
{
  delete ref;
  delete dv;
  delete sv;
}

// Copy from a std vector
void BlockVectorTest::testConstructor1()
{
  cout << "==================================" << endl;
  cout << "=== BlockVector tests start ...=== " << endl;
  cout << "==================================" << endl;
  cout << "--> Test: constructor 1." << endl;
  SiconosVector * w = new SimpleVector(3, 2);
  BlockVector * v = new BlockVector(*w); // Copy from a Simple

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->getNumberOfBlocks() == 1, true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", (*v)(i) == 2, true);
  Index tab = v->getTabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", tab.size() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", tab[0] == 3, true);
  SiconosVector * z = v->getVectorPtr(0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", z->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", z->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", z->getNum() == 1, true);
  delete v;
  delete w;
  cout << "--> Constructor 1 test ended with success." << endl;
}

void BlockVectorTest::testConstructor2()
{
  cout << "--> Test: constructor 2." << endl;
  SiconosVector * w = new SimpleVector(3, 2);
  SiconosVector * z = new SimpleVector(5, 3, SPARSE);
  BlockVector * x = new BlockVector(*w); // Copy from a Simple
  x->addPtr(z);
  BlockVector * v = new BlockVector(*x); // Copy from a Block
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v->size() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v->getNumberOfBlocks() == 2, true);
  for (unsigned int i = 0; i < w->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*v)(i) == 2, true);
  for (unsigned int i = w->size(); i < z->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*v)(i) == 3, true);
  Index tab = v->getTabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", tab.size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", tab[0] == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", tab[1] == 8, true);
  SiconosVector * ww = v->getVectorPtr(0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", ww->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", ww->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", ww->getNum() == 1, true);
  ww = v->getVectorPtr(1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", ww->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", ww->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", ww->getNum() == 4, true);
  delete v;
  delete w;
  delete x;
  delete z;
  cout << "--> Constructor 2 test ended with success." << endl;
}

// copy from a SiconosVector
void BlockVectorTest::testConstructor3()
{
  cout << "--> Test: constructor 3." << endl;
  SiconosVector * w = new SimpleVector(3, 2);
  SiconosVector * z = new SimpleVector(5, 3, SPARSE);
  BlockVector * v = new BlockVector(w, z);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->size() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->getNumberOfBlocks() == 2, true);
  for (unsigned int i = 0; i < w->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*v)(i) == 2, true);
  for (unsigned int i = w->size(); i < z->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*v)(i) == 3, true);
  Index tab = v->getTabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab.size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab[0] == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab[1] == 8, true);
  SiconosVector * ww = v->getVectorPtr(0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", ww == w, true);
  ww = v->getVectorPtr(1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", ww == z, true);
  delete v;
  delete w;
  delete z;
  cout << "--> Constructor 3 test ended with success." << endl;
}

// with number of blocks and their (common) size.
void BlockVectorTest::testConstructor4()
{
  cout << "--> Test: constructor 4." << endl;
  BlockVector * v = new BlockVector(3, 4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->size() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->getNumberOfBlocks() == 3, true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", (*v)(i) == 0, true);
  Index tab = v->getTabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", tab.size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", tab[0] == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", tab[1] == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", tab[2] == 12, true);
  delete v;
  cout << "--> Constructor 4 test ended with success." << endl;
}

// Constructor for block of blocks
void BlockVectorTest::testConstructor5()
{
  cout << "--> Test: constructor 5." << endl;
  SiconosVector * q = new BlockVector(3, 3);
  std::vector<SiconosVector*> in;
  in.push_back((*q)[0]);
  in.push_back((*q)[1]);
  in.push_back((*q)[1]);
  in.push_back((*q)[2]);
  BlockVector * v = new BlockVector(2, 2, in);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*v)[0]->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*v)[1]->isBlock(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->size() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*v)[0]->size() == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*v)[1]->size() == 6, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->getNumberOfBlocks() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*v)[0]->getNumberOfBlocks() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*v)[1]->getNumberOfBlocks() == 2, true);

  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*v)(i) == 0, true);
  Index tab = v->getTabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab.size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab[0] == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab[1] == 12, true);
  tab = (*v)[0]->getTabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab.size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab[0] == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab[1] == 6, true);
  tab = (*v)[1]->getTabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab.size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab[0] == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab[1] == 6, true);

  delete v;
  delete q;
  cout << "--> Constructor 5 test ended with success." << endl;
}

// zero
void BlockVectorTest::testZero()
{
  cout << "--> Test: zero." << endl;
  SiconosVector *v = new BlockVector(*ref);
  v->zero();
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*v)(i) == 0, true);
  delete v;
  cout << "--> zero test ended with success." << endl;
}

void BlockVectorTest::testNorm()
{
  cout << "--> Test: norm." << endl;
  SiconosVector * w = new SimpleVector(3, 2);
  SiconosVector * z = new SimpleVector(5, 3);
  BlockVector *v = new BlockVector(w, z);
  double n2 = v->norm2();
  v->addPtr(ref);
  double ni = v->normInf();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNorm : ", n2 == sqrt(57), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNorm : ", ni == ref->normInf(), true);
  delete v;
  delete w;
  delete z;
  cout << "--> norm test ended with success." << endl;
}

// OPERATORS

// =
void BlockVectorTest::testAssignment()
{
  cout << "--> Test: assignment." << endl;
  SiconosVector * w0 = new SimpleVector(3, 2);
  SiconosVector * z0 = new SimpleVector(5, 3);
  SiconosVector * w1 = new SimpleVector(3, 4);
  SiconosVector * z1 = new SimpleVector(5, 5);
  SiconosVector *v = new BlockVector(w0, z0);
  SiconosVector *x = new BlockVector(w1, z1);
  *v = *x; // assign a block to a block.
  SiconosVector * test = v->getVectorPtr(0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", *test == *w0, true);
  test = v->getVectorPtr(1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", *test == *z0, true);

  SiconosVector * w2 = new SimpleVector(8, 7);
  *v = *w2; // assign a simple to a block.
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == 7, true);

  delete w2;
  delete x;
  delete v;
  delete z1;
  delete w1;
  delete z0;
  delete w0;
  cout << "--> assignment test ended with success." << endl;
}

// +=
void BlockVectorTest::testOperators1()
{
  cout << "--> Test: operators1." << endl;
  SiconosVector *v = new BlockVector(2, 3);
  v->fill(4);
  SiconosVector *x = new BlockVector(2, 3);
  x->fill(5);
  *v += *x; // += a block.
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == 9, true);
  SiconosVector * z = new SimpleVector(6);
  z->fill(19);
  *v += *z;
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == 28, true);

  delete x;
  delete v;
  delete z;
  cout << "--> operators1 test ended with success." << endl;
}

// -=
void BlockVectorTest::testOperators2()
{
  cout << "--> Test: operators2." << endl;
  SiconosVector *v = new BlockVector(2, 3);
  v->fill(4);
  SiconosVector *x = new BlockVector(2, 3);
  x->fill(5);
  *v -= *x; // += a block.
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == -1, true);
  SiconosVector * z = new SimpleVector(6);
  z->fill(19);
  *v -= *z;
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == -20, true);

  delete x;
  delete v;
  delete z;
  cout << "--> operators2 test ended with success." << endl;
}

// *
void BlockVectorTest::testOperators3()
{
  cout << "--> Test: operators3." << endl;
  SiconosVector *v = new BlockVector(2, 3);
  v->fill(4);
  double multD = 2.3;
  int multI = 2;
  *v *= multD;
  *v *= multI;
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == multD * multI * 4, true);

  cout << "--> operators3 test ended with success." << endl;
}

// /
void BlockVectorTest::testOperators4()
{
  cout << "--> Test: operators4." << endl;
  SiconosVector *v = new BlockVector(2, 3);
  v->fill(4);
  double multD = 2.3;
  int multI = 2;
  *v /= multD;
  *v /= multI;
  double res = 4.0 / multD / multI;
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == res, true);
  cout << "--> operators4 test ended with success." << endl;
}

void BlockVectorTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of BlockVector Tests ===== " << endl;
  cout << "======================================" << endl;
}



