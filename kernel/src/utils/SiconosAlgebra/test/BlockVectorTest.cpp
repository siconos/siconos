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
#include "BlockVectorTest.hpp"

#include "BlockVector.hpp"

using namespace boost::numeric::ublas;
using BlockVector = siconos::algebra::BlockVector;

CPPUNIT_TEST_SUITE_REGISTRATION(BlockVectorTest);

void BlockVectorTest::setUp() {
  ref = std::make_shared<BlockVector>(1, 5);
  for (unsigned int i = 0; i < 5; ++i) (*ref)(i) = i;

  vq.resize(5, 1);
  for (unsigned int i = 0; i < 5; i++) vq[i] = i + 1;

  dv = std::make_shared<siconos::algebra::DenseVect>(3);
  (*dv)(0) = 1;
  (*dv)(1) = 2;
  (*dv)(2) = 3;
  sv = std::make_shared<siconos::algebra::SparseVect>(5);
  (*sv)(1) = 22;
}

void BlockVectorTest::tearDown() {}

// Copy from a std vector
void BlockVectorTest::testConstructor1() {
  std::cout << "==================================" << std::endl;
  std::cout << "=== BlockVector tests start ...=== " << std::endl;
  std::cout << "==================================" << std::endl;
  std::cout << "--> Test: constructor 1." << std::endl;
  auto w = std::make_shared<siconos::algebra::SiconosVector>(3, 2);
  // auto v = std::make_shared<BlockVector>(1, w->size());  // Copy from a Simple

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->size() == 3, true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->numberOfBlocks() == 1, true);
  // for (unsigned int i = 0; i < v->size(); i++)
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", (*v)(i) == 2, true);

  // auto tab = v->tabIndex();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", tab->size() == 1, true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", (*tab)[0] == 3, true);

  // auto z = v->vector(0);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", z->size() == 3, true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testConstructor1 : ", z->num() == siconos::algebra::UblasType::DENSE, true);
  // std::cout << "--> Constructor 1 test ended with success." << std::endl;
}

void BlockVectorTest::testConstructor2() {
  std::cout << "--> Test: constructor 2." << std::endl;
  auto w = std::make_shared<siconos::algebra::SiconosVector>(3, 2);
  // auto  z=
  // std::make_shared<siconos::algebra::SiconosVector>(5,3,siconos::algebra::UblasType::SPARSE);
  // Problem if z sparse: " Assertion failed in file
  // /usr/include/boost/numeric/ublas/vector_sparse.hpp at line 1253:
  // *this != (*this)().end() "

  auto z = std::make_shared<siconos::algebra::SiconosVector>(5, 3);
  auto x =
      std::make_shared<BlockVector>();  // Copy from a siconos::algebra::SiconosVector(Simple)
  x->insertPtr(w);
  x->insertPtr(z);

  auto v = std::make_shared<BlockVector>(*x);  // Copy from a Block

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v->size() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v->numberOfBlocks() == 2, true);
  for (unsigned int i = 0; i < w->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*v)(i) == 2, true);
  for (unsigned int i = w->size(); i < z->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*v)(i) == 3, true);
  auto tab = v->tabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", tab->size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*tab)[0] == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*tab)[1] == 8, true);
  auto ww = v->vector(0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", ww->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testConstructor2 : ", ww->num() == siconos::algebra::UblasType::DENSE, true);
  ww = v->vector(1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", ww->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testConstructor2 : ", ww->num() == siconos::algebra::UblasType::DENSE, true);

  auto x2 = std::make_shared<BlockVector>();
  x2->insertPtr(w);
  x2->insertPtr(z);
  auto v2 = std::make_shared<BlockVector>(*x2);  // Copy from a Siconos(Block)

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v2->size() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v2->numberOfBlocks() == 2, true);
  for (unsigned int i = 0; i < w->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*v2)(i) == 2, true);
  for (unsigned int i = w->size(); i < z->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*v2)(i) == 3, true);
  auto tab2 = v2->tabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", tab2->size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*tab2)[0] == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*tab2)[1] == 8, true);
  auto ww2 = v2->vector(0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", ww2->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testConstructor2 : ", ww2->num() == siconos::algebra::UblasType::DENSE, true);
  ww2 = v2->vector(1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", ww2->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testConstructor2 : ", ww2->num() == siconos::algebra::UblasType::DENSE, true);
  std::cout << "--> Constructor 2 test ended with success." << std::endl;
}

void BlockVectorTest::testConstructor3() {
  std::cout << "--> Test: constructor 3." << std::endl;
  auto w = std::make_shared<siconos::algebra::SiconosVector>(3, 2);
  auto z = std::make_shared<siconos::algebra::SiconosVector>(
      5, 3, siconos::algebra::UblasType::SPARSE);
  auto v = std::make_shared<BlockVector>(w, z);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->size() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->numberOfBlocks() == 2, true);
  for (unsigned int i = 0; i < w->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*v)(i) == 2, true);
  for (unsigned int i = w->size(); i < z->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*v)(i) == 3, true);
  auto tab = v->tabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab->size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[0] == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[1] == 8, true);
  auto ww = v->vector(0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", ww == w, true);
  ww = v->vector(1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", ww == z, true);
  std::cout << "--> Constructor 3 test ended with success." << std::endl;
}

// with number of blocks and their (common) size.
void BlockVectorTest::testConstructor4() {
  std::cout << "--> Test: constructor 4." << std::endl;
  auto v = std::make_shared<BlockVector>(3, 4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->size() == 12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->numberOfBlocks() == 3, true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", (*v)(i) == 0, true);
  auto tab = v->tabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", tab->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", (*tab)[0] == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", (*tab)[1] == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", (*tab)[2] == 12, true);
  std::cout << "--> Constructor 4 test ended with success." << std::endl;
}

// with number of blocks an empty vector
void BlockVectorTest::testConstructor5() {
  std::cout << "--> Test: constructor 5." << std::endl;
  auto v = std::make_shared<BlockVector>(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->size() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->numberOfBlocks() == 3, true);
  v->display();

  // test insertion
  auto w = std::make_shared<siconos::algebra::SiconosVector>(3, 2);
  v->insertPtr(w);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->numberOfBlocks() == 4, true);

  // set Vector
  v->setVectorPtr(2, w);
  v->display();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->size() == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->numberOfBlocks() == 4, true);

  auto ww = std::make_shared<siconos::algebra::SiconosVector>(5, 3);
  v->setVectorPtr(0, ww);
  v->display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->size() == 11, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->numberOfBlocks() == 4, true);

  v->insertPtr(ww);
  std::cout << "v->numberOfBlocks()" << v->numberOfBlocks() << std::endl;
  std::cout << "v->size()" << v->size() << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->size() == 16, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->numberOfBlocks() == 5, true);
  v->display();
  auto tab = v->tabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*tab)[0] == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*tab)[1] == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*tab)[2] == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*tab)[3] == 11, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*tab)[4] == 16, true);

  v->setVectorPtr(2, ww);
  v->display();
  tab = v->tabIndex();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->size() == 18, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->numberOfBlocks() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", tab->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*tab)[0] == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*tab)[1] == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*tab)[2] == 10, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*tab)[3] == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", (*tab)[4] == 18, true);

  std::cout << "--> Constructor 5 test ended with success." << std::endl;
}

// zero
void BlockVectorTest::testZero() {
  std::cout << "--> Test: zero." << std::endl;
  auto v = std::make_shared<BlockVector>(*ref);
  v->setZero();
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*v)(i) == 0, true);
  std::cout << "--> zero test ended with success." << std::endl;
}

void BlockVectorTest::testFill() {
  std::cout << "--> Test: fill." << std::endl;
  auto v = std::make_shared<BlockVector>(*ref);
  auto z = std::make_shared<siconos::algebra::SiconosVector>(5, 3);
  v->insertPtr(z);
  v->fill(4.3);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", (*v)(i) == 4.3, true);
  for (unsigned int i = 0; i < 5; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", (*z)(i) == 4.3, true);
  std::cout << "--> fill test ended with success." << std::endl;
}

void BlockVectorTest::testNorm() {
  //  std::cout << "--> Test: norm." <<std::endl;
  //  auto w= std::make_shared<siconos::algebra::SiconosVector>(3);
  //  auto z= std::make_shared<siconos::algebra::SiconosVector>(5);
  //
  //  (*w)(0)=1;
  //  (*w)(1)=2;
  //  (*w)(2)=3;
  //
  //  (*z)(0)=1;
  //  (*z)(1)=2;
  //  (*z)(2)=3;
  //  (*z)(3)=4;
  //  (*z)(4)=5;
  //
  //
  //  auto v= std::make_shared<BlockVector>(w,z);
  //  double n2 = v->norm2();
  //  v->insertPtr(ref.vector(0));
  ////  double ni = v->normInf();
  //
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNorm : ",
  //  n2-sqrt(69)<std::numeric_limits<double>::epsilon(), true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNorm : ", ni==5, true);
  std::cout << "--> norm test ended with success." << std::endl;
}

void BlockVectorTest::testSetBlock() {
  std::cout << "--> Test: set block in BlockVector." << std::endl;

  // Block copy from a Simple into a Block.
  auto ref2 = std::make_shared<siconos::algebra::SiconosVector>(3);
  for (unsigned int i = 0; i < 3; ++i) (*ref2)(i) = i;
  auto ref3 = std::make_shared<siconos::algebra::SiconosVector>(6);
  for (unsigned int i = 0; i < 6; ++i) (*ref3)(i) = i;
  auto BSV = std::make_shared<BlockVector>();
  BSV->insertPtr(ref3);
  BSV->insertPtr(ref2);

  auto vIn = std::make_shared<siconos::algebra::SiconosVector>(8);

  unsigned int sizeB = 6;
  unsigned int posIn = 1;
  unsigned int posOut = 2;

  BSV->setBlock(*vIn, sizeB, posIn, posOut);
  for (unsigned int i = 0; i < sizeB; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test setBlock : ", (*BSV)(i + posOut) == (*vIn)(i + posIn),
                                 true);

  std::cout << "--> setBlock3 test ended with success." << std::endl;
}
// OPERATORS

// =
void BlockVectorTest::testAssignment() {
  std::cout << "--> Test: assignment." << std::endl;
  auto w0 = std::make_shared<siconos::algebra::SiconosVector>(3, 2);
  auto z0 = std::make_shared<siconos::algebra::SiconosVector>(5, 3);
  auto w1 = std::make_shared<siconos::algebra::SiconosVector>(3, 4);
  auto z1 = std::make_shared<siconos::algebra::SiconosVector>(5, 5);

  auto v = std::make_shared<BlockVector>(w0, z0);

  // Block = Block (same sub-blocks sizes)
  auto xB = std::make_shared<BlockVector>(w1, z1);
  *v = *xB;
  auto test = v->vector(0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", *test == *w1, true);
  test = v->vector(1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", *test == *z1, true);
  // Block = Block (different sub-blocks sizes)
  xB = std::make_shared<BlockVector>(z1, w1);
  *v = *xB;
  for (unsigned int i = 0; i < v->size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == (*xB)(i), true);

  // Block = Siconos(Block) (same sub-blocks sizes)
  auto x = std::make_shared<BlockVector>(w1, z1);
  *v = *x;
  test = v->vector(0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", *test == *w1, true);
  test = v->vector(1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", *test == *z1, true);
  // Block = Siconos(Block) (different sub-blocks sizes)
  x = std::make_shared<BlockVector>(z1, w1);
  *v = *x;
  for (unsigned int i = 0; i < v->size(); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == (*x)(i), true);

  // Block = Siconos(Simple)
  auto w2 = std::make_shared<siconos::algebra::SiconosVector>(8, 7);

  *v = *w2;  // assign a simple to a block.
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == 7, true);
  std::cout << "--> assignment test ended with success." << std::endl;
}

// +=
void BlockVectorTest::testOperators1() {
  std::cout << "--> Test: operators1." << std::endl;
  unsigned int size1 = 3;
  unsigned int size2 = 2;
  auto x = std::make_shared<siconos::algebra::SiconosVector>(size1 + size2);
  x->fill(3.2);
  auto tmp1 = std::make_shared<siconos::algebra::SiconosVector>(size1);
  auto tmp2 = std::make_shared<siconos::algebra::SiconosVector>(size2);
  for (unsigned int i = 0; i < size1; ++i) (*tmp1)(i) = i;
  for (unsigned int i = 0; i < size2; ++i) (*tmp2)(i) = 100 * i;

  auto xB = std::make_shared<BlockVector>(tmp1, tmp2);

  auto v = std::make_shared<BlockVector>(*xB);  // Copy
  v->fill(4);
  // Block += Simple
  //  *v += *x;
  //  for (unsigned int i=0; i< size1+size2; i++)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE(" testOperators1: ", fabs((*v)(i) - 4 -(*x)(i))<tol,
  //    true);

  v->fill(4);
  // Block += Block
  *v += *xB;
  for (unsigned int i = 0; i < size1 + size2; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE(" testOperators1: ", fabs((*v)(i)-4 - (*xB)(i)) < tol, true);
  std::cout << "--> operators1 test ended with success." << std::endl;
}

// -=
void BlockVectorTest::testOperators2() {
  std::cout << "--> Test: operators2." << std::endl;
  unsigned int size1 = 3;
  unsigned int size2 = 2;
  auto x = std::make_shared<siconos::algebra::SiconosVector>(size1 + size2);
  x->fill(3.2);
  auto tmp1 = std::make_shared<siconos::algebra::SiconosVector>(size1);
  auto tmp2 = std::make_shared<siconos::algebra::SiconosVector>(size2);
  for (unsigned int i = 0; i < size1; ++i) (*tmp1)(i) = i;
  for (unsigned int i = 0; i < size2; ++i) (*tmp2)(i) = 100 * i;

  auto xB = std::make_shared<BlockVector>(tmp1, tmp2);

  auto v = std::make_shared<BlockVector>(*xB);  // Copy
  v->fill(4);
  // Block += Simple
  //  *v -= *x;
  //  for (unsigned int i=0; i< size1+size2; i++)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE(" testOperators2: ", fabs((*v)(i) - 4 +(*x)(i))<tol,
  //    true);

  v->fill(4);
  // Block += Block
  *v -= *xB;
  for (unsigned int i = 0; i < size1 + size2; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE(" testOperators2: ", fabs((*v)(i)-4 + (*xB)(i)) < tol, true);

  std::cout << "--> operators2 test ended with success." << std::endl;
}

// *
void BlockVectorTest::testOperators3() {
  std::cout << "--> Test: operators3." << std::endl;
  auto v = std::make_shared<BlockVector>(2, 3);
  v->fill(4);
  double multD = 2.3;
  int multI = 2;
  *v *= multD;
  *v *= multI;
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", (*v)(i) == multD * multI * 4, true);

  std::cout << "--> operators3 test ended with success." << std::endl;
}

// /=
void BlockVectorTest::testOperators4() {
  std::cout << "--> Test: operators4." << std::endl;
  auto v = std::make_shared<BlockVector>(2, 3);
  v->fill(4);
  double multD = 2.3;
  int multI = 2;
  *v /= multD;
  *v /= multI;
  double res = 4.0 / multD / multI;
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE(" testOperators4: ", (*v)(i) == res, true);
  std::cout << "--> operators4 test ended with success." << std::endl;
}

void BlockVectorTest::testInsert() {
  std::cout << "--> Test: insert." << std::endl;
  unsigned int size1 = 3, size2 = 2, size = size1 + size2;
  auto tmp1 = std::make_shared<siconos::algebra::SiconosVector>(size1);
  auto tmp2 = std::make_shared<siconos::algebra::SiconosVector>(size2);
  for (unsigned int i = 0; i < size1; ++i) (*tmp1)(i) = i;
  for (unsigned int i = 0; i < size2; ++i) (*tmp2)(i) = 100 * i;

  auto xB = std::make_shared<BlockVector>(tmp1, tmp2);

  auto y = std::make_shared<siconos::algebra::SiconosVector>(7);

  xB->insertPtr(y);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInsert : ", xB->numberOfBlocks() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInsert : ", xB->size() == (7 + size), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInsert : ", xB->vector(2) == y, true);

  xB->insertPtr(y);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInsert : ", xB->numberOfBlocks() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInsert : ", xB->size() == (14 + size), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInsert : ", xB->vector(3) == y, true);

  std::cout << "--> insert test ended with success." << std::endl;
}

void BlockVectorTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of BlockVector Tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
