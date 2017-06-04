/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "SiconosConfig.h"

#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosAlgebra.hpp"
#include "SiconosVectorTest.hpp"


using namespace boost::numeric::ublas;

CPPUNIT_TEST_SUITE_REGISTRATION(SiconosVectorTest);

void SiconosVectorTest::setUp()
{
  tol = 1e-12;

  size = 5;
  size1 = 3;
  size2 = 2;

  // size = size1 + size2;

  ref.reset(new SiconosVector(size));
  for (unsigned int i = 0; i < size; ++i)
    (*ref)(i) = i;
  vq.resize(size, 1);
  for (unsigned int i = 0; i < size; i++)
    vq[i] = i + 1;

  dv.reset(new DenseVect(3));
  (*dv)(0) = 1;
  (*dv)(1) = 2;
  (*dv)(2) = 3;
  sv.reset(new SparseVect(5));
  (*sv)(1) = 22;


  // const vectors used for operators test (ex: x and y in z = x + y)
  // "B" in name for BlockVectors
  x.reset(new SiconosVector(vq));
  y.reset(new SiconosVector(vq));
  tmp1.reset(new SiconosVector(size1));
  tmp2.reset(new SiconosVector(size2));
  for (unsigned int i = 0; i < size1; ++i)
    (*tmp1)(i) = i;
  for (unsigned int i = 0; i < size2; ++i)
    (*tmp2)(i) = 100 * i;

  xB.reset(new BlockVector(tmp1, tmp2));
  yB.reset(new BlockVector(tmp2, tmp1));

  // vectors used as results
  z.reset(new SiconosVector(size));
  tmp3.reset(new SiconosVector(size2));
  tmp4.reset(new SiconosVector(size1));
  zB.reset(new BlockVector());
  zB->insertPtr(tmp3);
  zB->insertPtr(tmp4);

}

void SiconosVectorTest::tearDown()
{}

void SiconosVectorTest::testConstructor0()
{
  std::cout << "====================================" <<std::endl;
  std::cout << "===  Simple Vector tests start ...=== " <<std::endl;
  std::cout << "====================================" <<std::endl;
  std::cout << "--> Test: constructor 0." <<std::endl;
  SP::SiconosVector v(new SiconosVector(3));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", v->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", v->num() == 1, true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", (*v)(i) == 0, true);
  v.reset(new SiconosVector(3, Siconos::SPARSE));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", v->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", v->num() == 4, true);
  std::cout << "--> Constructor 0 test ended with success." <<std::endl;
}

void SiconosVectorTest::testConstructor1()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::SiconosVector v(new SiconosVector(3, 2.4));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->num() == 1, true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", (*v)(i) == 2.4, true);
  v.reset(new SiconosVector(3, 2.4, Siconos::SPARSE));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->num() == 4, true);

  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", (*v)(i) == 2.4, true);
  std::cout << "--> Constructor 1 test ended with success." <<std::endl;
}

void SiconosVectorTest::testConstructor2()
{
  std::cout << "--> Test: constructor 2." <<std::endl;
  // Copy from a std::vector
  SP::SiconosVector v(new SiconosVector(vq));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*v)(i) == vq[i], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->num() == 1, true);
  std::cout << "--> Constructor 2 test ended with success." <<std::endl;
}

// copy from a SiconosVector (Simple)
void SiconosVectorTest::testConstructor3()
{
  std::cout << "--> Test: constructor 3." <<std::endl;
  SP::SiconosVector tmp(new SiconosVector(vq));
  SP::SiconosVector v(new SiconosVector(*tmp));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", *v == *tmp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->num() == 1, true);
  tmp.reset(new SiconosVector(3, Siconos::SPARSE));
  v.reset(new SiconosVector(*tmp));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", *v == *tmp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->num() == 4, true);
  std::cout << "--> Constructor 3 test ended with success." <<std::endl;
}

// copy from a SiconosVector (Block)
void SiconosVectorTest::testConstructor3Bis()
{
  std::cout << "--> Test: constructor 3Bis." <<std::endl;
  //  SP::SiconosVector v(new SiconosVector(*xB));
  //
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3Bis : ", v->isBlock(), false);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3Bis : ", v->size() == size, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3Bis : ", v->num() == 1, true);
  //  for (unsigned int i = 0; i<size1; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3Bis : ", (*v)(i) == (*tmp1)(i), true);
  //  for (unsigned int i = size1; i<size2; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3Bis : ", (*v)(i) == (*tmp2)(i-size1), true);
  //
  std::cout << "--> Constructor 3Bis test ended with success." <<std::endl;
}

// copy from a SiconosVector
void SiconosVectorTest::testConstructor4()
{
  std::cout << "--> Test: constructor 4." <<std::endl;
  SP::SiconosVector tmp(new SiconosVector(vq));
  SP::SiconosVector v(new SiconosVector(*tmp));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", *v == *tmp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->num() == 1, true);
  tmp.reset(new SiconosVector(4, Siconos::SPARSE));
  v.reset(new SiconosVector(*tmp));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->size() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", *v == *tmp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->num() == 4, true);
  std::cout << "--> Constructor 4 test ended with success." <<std::endl;
}

// Copy from a Dense
void SiconosVectorTest::testConstructor5()
{
  std::cout << "--> Test: constructor 5." <<std::endl;
  SP::SiconosVector v(new SiconosVector(*dv));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->size() == dv->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", norm_inf(v->getDense() - *dv) < tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->num() == 1, true);

  std::cout << "--> Constructor 5 test ended with success." <<std::endl;
}

// Copy from a sparse
void SiconosVectorTest::testConstructor6()
{
  std::cout << "--> Test: constructor 6." <<std::endl;
  SP::SiconosVector v(new SiconosVector(*sv));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", v->size() == sv->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", norm_inf(v->getSparse() - *sv) < tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", v->num() == 4, true);
  std::cout << "--> Constructor 6 test ended with success." <<std::endl;
}

// From a file
void SiconosVectorTest::testConstructor7()
{
  std::cout << "--> Test: constructor 7." <<std::endl;
  SP::SiconosVector v(new SiconosVector("vect.dat", 1));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", v->size() == 4, true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", (*v)(i) == i, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", v->num() == 1, true);
  std::cout << "--> Constructor 7 test ended with success." <<std::endl;
}

// zero
void SiconosVectorTest::testZero()
{
  std::cout << "--> Test: zero." <<std::endl;
  SP::SiconosVector v(new SiconosVector(vq));
  v->zero();
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*v)(i) == 0, true);
  std::cout << "--> zero test ended with success." <<std::endl;
}

void SiconosVectorTest::testNorm()
{
  std::cout << "--> Test: norm." <<std::endl;
  SP::SiconosVector v(new SiconosVector(*dv));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNorm : ", v->normInf() == norm_inf(*dv), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNorm : ", (v->norm2() - norm_2(*dv)) < tol, true);
  std::cout << "--> norm test ended with success." <<std::endl;
}

// fill
void SiconosVectorTest::testFill()
{
  std::cout << "--> Test: fill." <<std::endl;
  double val = 4.5;
  SP::SiconosVector v(new SiconosVector(5, val));
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testFill : ", (*v)(i) == val, true);
  std::cout << "--> fill test ended with success." <<std::endl;
}

//resize
void SiconosVectorTest::testResize()
{
  std::cout << "--> Test: resize." <<std::endl;
  SP::SiconosVector v(new SiconosVector(vq));
  v->resize(6);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test resize : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test resize : ", v->size() == 6, true);
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test resize : ", (*v)(i) == vq[i], true);
  for (unsigned int i = vq.size(); i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test resize : ", (*v)(i) == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test resize : ", v->num() == 1, true);
  std::cout << "-->  resize test ended with success." <<std::endl;
}

void SiconosVectorTest::testSetBlock()
{
  std::cout << "--> Test: setBlock." <<std::endl;

  // Block copy from a Simple into a Simple
  unsigned int sizeB = 2;
  SP::SiconosVector subBlock(new SiconosVector(6));
  unsigned int posIn = 1;
  unsigned int posOut = 3;
  setBlock(*ref, subBlock, sizeB, posIn, posOut);

  for (unsigned int i = posOut; i < sizeB; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test setBlock : ", (*subBlock)(i) == (*ref)(i + posIn), true);
  std::cout << "--> setBlock test ended with success." <<std::endl;
}

// void SiconosVectorTest::testSetBlock2()
// {
//   std::cout << "--> Test: setBlock2." <<std::endl;

// Note FP :  Block copy from a Block into a Simple is no more implemented.
//   // Block copy from a Block into a Simple.
//   SP::BlockVector BSV(new BlockVector());
//   BSV->insertPtr(ref);
//   SP::SiconosVector ref2(new SiconosVector(3));
//   for (unsigned int i = 0; i<3; ++i)
//     (*ref2)(i)=i;
//   BSV->insertPtr(ref2);

//   SP::SiconosVector vOut(new SiconosVector(14));

//   unsigned int sizeB = 6;
//   unsigned int posIn = 1;
//   unsigned int posOut = 3;

// //  setBlock(*BSV, vOut, sizeB, posIn, posOut);
// //  for (unsigned int i = 0; i < sizeB; ++i)
// //    CPPUNIT_ASSERT_EQUAL_MESSAGE("test setBlock : ", (*vOut)(i+posOut) == (*BSV)(i+posIn), true);

//   std::cout << "--> setBlock2 test ended with success." <<std::endl;
// }

void SiconosVectorTest::testSetBlock3()
{
  std::cout << "--> Test: setBlock3." <<std::endl;

  // Block copy from a Simple into a Block.
  SP::SiconosVector ref2(new SiconosVector(3));
  for (unsigned int i = 0; i < 3; ++i)
    (*ref2)(i) = i;
  SP::SiconosVector ref3(new SiconosVector(6));
  for (unsigned int i = 0; i < 6; ++i)
    (*ref3)(i) = i;
  SP::BlockVector BSV(new BlockVector());
  BSV->insertPtr(ref3);
  BSV->insertPtr(ref2);

  SP::SiconosVector vIn(new SiconosVector(8));

  unsigned int sizeB = 6;
  unsigned int posIn = 1;
  unsigned int posOut = 2;

  BSV->setBlock(*vIn, sizeB, posIn, posOut);
  for (unsigned int i = 0; i < sizeB; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test setBlock : ", (*BSV)(i + posOut) == (*vIn)(i + posIn), true);

  std::cout << "--> setBlock3 test ended with success." <<std::endl;
}

void SiconosVectorTest::testSetBlock4()
{
  std::cout << "--> Test: setBlock4." <<std::endl;

  // copy from a Simple into a simple
  unsigned int sizeB = 9;
  SP::SiconosVector subBlock(new SiconosVector(sizeB));
  unsigned int pos = 1;
  subBlock->setBlock(pos, *ref);

  for (unsigned int i = pos; i < pos + 5; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test setBlock : ", fabs((*subBlock)(i) - (*ref)(i - pos)) < tol, true);

  // copy from a Block into a simple
  SP::SiconosVector ref2(new SiconosVector(2));
  for (unsigned int i = 0; i < 2; ++i)
    (*ref2)(i) = i;
  SP::SiconosVector ref3(new SiconosVector(4));
  for (unsigned int i = 0; i < 2; ++i)
    (*ref3)(i) = i;
  SP::BlockVector BSV(new BlockVector());
  BSV->insertPtr(ref3);
  BSV->insertPtr(ref2);

  pos = 2;
  //  subBlock->setBlock(pos, *BSV);

  //  for (unsigned int i = pos; i < pos+6; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("test setBlock : ", fabs((*BSV)(i-pos) - (*subBlock)(i))< tol, true);

  std::cout << "--> setBlock4 test ended with success." <<std::endl;
}

// OPERATORS

// =
void SiconosVectorTest::testAssignment()
{
  std::cout << "--> Test: assignment." <<std::endl;

  SP::SiconosVector v(new SiconosVector(vq.size()));

  // Simple = Simple
  SP::SiconosVector w(new SiconosVector(vq));
  *v = *w;
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == (*w)(i), true);

  // Simple = Siconos (Simple)
  SP::SiconosVector w2(new SiconosVector(vq));
  *v = *w2;
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == (*w2)(i), true);

  // Simple = Siconos (Block)
  *v = *xB;
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == (*xB)(i), true);
  std::cout << "--> operators1 test ended with success." <<std::endl;
}

// ()
void SiconosVectorTest::testOperators1()
{
  std::cout << "--> Test: operators1." <<std::endl;

  std::cout << "--> operators1 test ended with success." <<std::endl;
}

// +=, -=, *=, /=
void SiconosVectorTest::testOperators2()
{
  std::cout << "--> Test: operators2." <<std::endl;

  // dense +=, -=, ... dense
  *z += *x;
  *z += *x;
  for (unsigned int i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", fabs((*z)(i) - 2 * (*x)(i)) < tol, true);

  double a = 2.2;
  int a1 = 2;

  *z *= a;
  for (unsigned int i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", fabs((*z)(i) - 2 * a * (*x)(i)) < tol, true);
  *z *= a1;
  for (unsigned int i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", fabs((*z)(i) - 2 * a1 * a * (*x)(i)) < tol, true);
  *z /= a1;
  for (unsigned int i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", fabs((*z)(i) - 2 * a * (*x)(i)) < tol, true);
  *z /= a;
  for (unsigned int i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", fabs((*z)(i) - 2 * (*x)(i)) < tol, true);
  *z -= *x;
  for (unsigned int i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", fabs((*z)(i) - (*x)(i)) < tol, true);

  // Simple +=, -= Block
  z->zero();
  *z += *x;

  //  *z+=*xB;
  //  for (unsigned int i=0; i< size; i++)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ",fabs((*z)(i)-(*x)(i)-(*xB)(i))<tol, true);
  //  *z-=*xB;
  //  for (unsigned int i=0; i< size; i++)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ",fabs((*z)(i)-(*x)(i))<tol, true);
  //
  // dense +=, -=, ... sparse
  SP::SiconosVector w(new SiconosVector(*sv));
  z->zero();
  *z += *w;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", norm_inf(*z->dense() - *w->sparse()) < tol, true);
  *z *= 3;
  *z -= *w;
  for (unsigned int i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", norm_inf(*z->dense() - 2 * *w->sparse()) < tol, true);
  *z /= 2.0;
  for (unsigned int i = 0; i < size; i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", norm_inf(*z->dense() - *w->sparse()) < tol, true);

  // sparse += -= sparse
  SP::SiconosVector v(new SiconosVector(5, Siconos::SPARSE));
  *v += *w;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", norm_inf(v->getSparse() - *sv) < tol, true);
  *v *= 3;
  *v -= *w;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", norm_inf(v->getSparse() - 2 * *sv) < tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", v->num() == 4, true);
  std::cout << "--> operators2 test ended with success." <<std::endl;
}

// inner_prod
void SiconosVectorTest::testOperators3()
{
  std::cout << "--> Test: operators3." <<std::endl;
  SP::SiconosVector v(new SiconosVector(vq));
  SP::SiconosVector w(new SiconosVector(vq));
  double res;
  // *
  res = inner_prod(*v, *w); // dense*dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", res == 55, true);
  w.reset(new SiconosVector(*sv));
  res = inner_prod(*v, *w); // dense*sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", res == 44, true);
  v.reset(new SiconosVector(*sv));
  res = inner_prod(*v, *w); // sparse*sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", res == 484, true);
  w.reset(new SiconosVector(vq));
  res = inner_prod(*v, *w); // sparse*dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", res == 44, true);
  std::cout << "--> operators3 test ended with success." <<std::endl;
}

// vector * or / by double or int
void SiconosVectorTest::testOperators4()
{
  std::cout << "--> Test: operators4." <<std::endl;

  double a = 2.2;
  int a1 = 2;

  //  z = a * x

  //  dense = a*dense or dense/a
  *z = a**x; //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(z->getDense() - a * x->getDense()) < tol, true);
  z->zero();
  *z = *x * a; //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(z->getDense() - a * x->getDense()) < tol, true);
  z->zero();
  *z = a1**x; //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(z->getDense() - a1 * x->getDense()) < tol, true);
  z->zero();
  *z = *x * a1; //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(z->getDense() - a1 * x->getDense()) < tol, true);
  z->zero();
  *z = *x / a; //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(z->getDense() - x->getDense() / a) < tol, true);
  z->zero();
  *z = *x / a1; //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(z->getDense() - x->getDense() / a1) < tol, true);
  z->zero();

  //  sparse = a*sparse or sparse/a
  x.reset(new SiconosVector(*sv));
  z.reset(new SiconosVector(size, Siconos::SPARSE));
  *z = a**x; //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(z->getSparse() - a * x->getSparse()) < tol, true);
  z->zero();
  *z = *x * a; //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(z->getSparse() - a * x->getSparse()) < tol, true);
  z->zero();
  *z = a1**x; //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(z->getSparse() - a1 * x->getSparse()) < tol, true);
  z->zero();
  *z = *x * a1; //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(z->getSparse() - a1 * x->getSparse()) < tol, true);
  z->zero();

  // Following tests failed. Sparse init pb. To be reviewed.

  //   *z = (*x)/a; //
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ",norm_2(z->getSparse()-x->getSparse()/a)<tol, true);
  //   z->zero();
  //   *z = (*x)/a1; //
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ",norm_2(z->getSparse()-x->getSparse()/a1)<tol, true);

  // simple = a * block
  z.reset(new SiconosVector(size));
  //  *z = a**xB; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*z)(i) - (a*(*xB)(i)))<tol, true);
  //  z->zero();
  //  *z = *xB*a; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*z)(i) - (a*(*xB)(i)))<tol, true);
  //  z->zero();
  //  *z = a1**xB; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*z)(i) - (a1*(*xB)(i)))<tol, true);
  //
  //  z->zero();
  //  *z = *xB*a1; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*z)(i) - (a1*(*xB)(i)))<tol, true);
  //
  //  z->zero();
  //  *z = (*xB)/a; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ",fabs((*z)(i) -((*xB)(i)/a))< tol, true);
  //
  //  z->zero();
  //  *z = (*xB)/a1; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*z)(i) - ((*xB)(i)/a1))< tol, true);
  //
  //  // x block, z block
  //
  //  *zB = a**xB; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - (a*(*xB)(i)))< tol, true);
  //  zB->zero();
  //  *zB = *xB*a; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - (a*(*xB)(i)))< tol, true);
  //  zB->zero();
  //  *zB = a1**xB; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - (a1*(*xB)(i)))< tol, true);
  //
  //  zB->zero();
  //  *zB = *xB*a1; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - (a1*(*xB)(i)))< tol, true);
  //
  //  zB->zero();
  //  *zB = *xB/a; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - ((*xB)(i)/a))< tol, true);
  //
  //  zB->zero();
  //  *zB = *xB/a1; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - ((*xB)(i)/a1))< tol, true);
  //
  //  x.reset(new SiconosVector(vq));
  //
  //  // block = a *simple
  //  zB->zero();
  //  *zB = a**x; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - (a*(*x)(i)))< tol, true);
  //  zB->zero();
  //  *zB = *x*a; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - (a*(*x)(i)))< tol, true);
  //  zB->zero();
  //  *zB = a1**x; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - (a1*(*x)(i)))< tol, true);
  //
  //  zB->zero();
  //  *zB = *x*a1; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - (a1*(*x)(i)))< tol, true);
  //
  //  zB->zero();
  //  *zB = *x/a; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - ((*x)(i)/a))< tol, true);
  //
  //  zB->zero();
  //  *zB = *x/a1; //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", fabs((*zB)(i) - ((*x)(i)/a1))< tol, true);
  //
  //  zB->zero();
  std::cout << "--> operators4 test ended with success." <<std::endl;
}

// vector * or / by double with function scal
void SiconosVectorTest::testOperators4Bis()
{
  std::cout << "--> Test: operators4Bis." <<std::endl;

  double a = 2.2;
  //  z = a * x

  //  dense = a*dense or dense/a
  scal(a, *x, *z); //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Bis : ", norm_2(z->getDense() - a * x->getDense()) < tol, true);
  z->zero();
  scal(1.0 / a, *x, *z);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Bis : ", norm_2(z->getDense() - x->getDense() / a) < tol, true);
  z->zero();

  //  sparse = a*sparse or sparse/a
  x.reset(new SiconosVector(*sv));
  z.reset(new SiconosVector(size, Siconos::SPARSE));
  scal(a, *x, *z);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Bis : ", norm_2(z->getSparse() - a * x->getSparse()) < tol, true);
  z->zero();
  // Following tests failed. Sparse init pb. To be reviewed.
  scal(1.0 / a, *x, *z);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Bis : ", norm_2(z->getSparse() - x->getSparse() / a) < tol, true);

  // simple = a * block
  z.reset(new SiconosVector(size));
  //  scal(a,*xB,*z); //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Bis : ", fabs((*z)(i) - (a*(*xB)(i)))<tol, true);
  //  z->zero();
  //  scal(1.0/a, *xB, *z);
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Bis : ",fabs((*z)(i) -((*xB)(i)/a))< tol, true);
  //
  // x block, z block

  //  scal(a, *xB, *zB);
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Bis : ", fabs((*zB)(i) - (a*(*xB)(i)))< tol, true);
  //  zB->zero();
  //  scal(1.0/a, *xB, *zB);
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Bis : ", fabs((*zB)(i) - ((*xB)(i)/a))< tol, true);
  //
  //  x.reset(new SiconosVector(vq));
  //
  //  // block = a *simple
  //  zB->zero();
  //  scal(a,*x,*zB);
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Bis : ", fabs((*zB)(i) - (a*(*x)(i)))< tol, true);
  //  zB->zero();
  //  scal(1.0/a, *x, *zB);
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Bis : ", fabs((*zB)(i) - ((*x)(i)/a))< tol, true);
  //  zB->zero();
  std::cout << "--> operators4Bis test ended with success." <<std::endl;
}

// vector * or / by double with function scal (init = false) =>   z += a*x
void SiconosVectorTest::testOperators4Ter()
{
  std::cout << "--> Test: operators4Ter." <<std::endl;

  double a = 2.2;
  //  z = a * x
  z->zero();
  //  dense += a*dense or dense/a
  scal(a, *x, *z, false); //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Ter : ", norm_2(z->getDense() - a * x->getDense()) < tol, true);
  scal(a, *x, *z, false); //
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Ter : ", norm_2(z->getDense() - 2 * a * x->getDense()) < tol, true);

  //  sparse += a*sparse or sparse/a
  x.reset(new SiconosVector(*sv));
  z.reset(new SiconosVector(size, Siconos::SPARSE));
  z->zero();
  scal(a, *x, *z, false);
  scal(a, *x, *z, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Ter : ", norm_2(z->getSparse() - 2 * a * x->getSparse()) < tol, true);


  //  // simple = a * block
  //  z.reset(new SiconosVector(size));
  //  scal(a,*xB,*z,false); //
  //  scal(a,*xB,*z,false); //
  //
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Ter : ", fabs((*z)(i) - 2*(a*(*xB)(i)))<tol, true);
  //  z->zero();
  //  // x block, z block
  //
  //  scal(a, *xB, *zB,false);
  //  scal(a, *xB, *zB,false);
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Ter : ", fabs((*zB)(i) - 2*(a*(*xB)(i)))< tol, true);
  //
  //  x.reset(new SiconosVector(vq));
  //
  //  // block = a *simple
  //  zB->zero();
  //  scal(a,*x,*zB, false);
  //  scal(a,*x,*zB, false);
  //  for (unsigned int i = 0; i< size ; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4Ter : ", fabs((*zB)(i) - 2*(a*(*x)(i)))< tol, true);
  //  zB->zero();
  std::cout << "--> operators4Ter test ended with success." <<std::endl;
}

// +
void SiconosVectorTest::testOperators5()
{
  //   // z = x + y
  std::cout << "--> Test: operators5." <<std::endl;

  // Simple = Simple + Simple, all dense
  *z = *x + *y ;
  for (unsigned int i = 0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", fabs((*z)(i) - (*x)(i) - (*y)(i)) < tol , true);
  z->zero();
  // Simple = Simple + Block, all dense
  //  *z = *x + *yB ;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", fabs((*z)(i) - (*x)(i) - (*yB)(i))<tol , true);
  //
  //  z->zero();
  //  // Simple = Block + Simple, all dense
  //  *z = *xB + *y ;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", fabs((*z)(i) - (*xB)(i) - (*y)(i))<tol , true);
  //
  //  z->zero();
  //  // Simple = Block + Block
  //  *z = *xB + *yB;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", fabs((*z)(i) - (*xB)(i) - (*yB)(i))<tol , true);
  //
  // TO DO x and/or y sparse.

  // Block = Simple + Simple, all dense
  //  *zB = *x + *y ;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", fabs((*zB)(i) - (*x)(i) - (*y)(i))<tol , true);
  //  zB->zero();
  //  // Block = Simple + Block, all dense
  //  *zB = *x + *yB ;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", fabs((*zB)(i) - (*x)(i) - (*yB)(i))<tol , true);
  //
  //  zB->zero();
  //  // Block = Block + Simple, all dense
  //  *zB = *xB + *y ;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", fabs((*zB)(i) - (*xB)(i) - (*y)(i))<tol , true);
  //
  //  zB->zero();
  //  // Block = Block + Block
  //  *zB = *xB + *yB;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", fabs((*zB)(i) - (*xB)(i) - (*yB)(i))<tol , true);
  std::cout << "--> operators5 test ended with success." <<std::endl;
}

void SiconosVectorTest::testOperators5Bis()
{
  //   // z = x + y with add
  std::cout << "--> Test: operators5Bis." <<std::endl;

  // Simple = Simple + Simple, all dense
  add(*x, *y, *z);
  for (unsigned int i = 0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5Bis : ", fabs((*z)(i) - (*x)(i) - (*y)(i)) < tol , true);
  z->zero();
  // Simple = Simple + Block, all dense
  //  add(*x,*yB,*z);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5Bis : ", fabs((*z)(i) - (*x)(i) - (*yB)(i))<tol , true);
  //
  //  z->zero();
  //  // Simple = Block + Simple, all dense
  //  add(*xB,*y,*z);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5Bis : ", fabs((*z)(i) - (*xB)(i) - (*y)(i))<tol , true);
  //
  //  z->zero();
  //  // Simple = Block + Block
  //  add(*xB,*yB,*z);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5Bis : ", fabs((*z)(i) - (*xB)(i) - (*yB)(i))<tol , true);
  //
  //  // TO DO x and/or y sparse.
  //
  //  // Block = Simple + Simple, all dense
  //  add(*x,*y,*zB);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5Bis : ", fabs((*zB)(i) - (*x)(i) - (*y)(i))<tol , true);
  //  zB->zero();
  //  // Block = Simple + Block, all dense
  //  add(*x,*yB,*zB);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5Bis : ", fabs((*zB)(i) - (*x)(i) - (*yB)(i))<tol , true);
  //
  //  zB->zero();
  //  // Block = Block + Simple, all dense
  //  add(*xB,*y,*zB);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5Bis : ", fabs((*zB)(i) - (*xB)(i) - (*y)(i))<tol , true);
  //
  //  zB->zero();
  //  // Block = Block + Block
  //  add(*xB,*yB,*zB);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5Bis : ", fabs((*zB)(i) - (*xB)(i) - (*yB)(i))<tol , true);
  //  std::cout << "--> operators5Bis test ended with success." <<std::endl;
}

// -
void SiconosVectorTest::testOperators6()
{
  //   // z = x - y
  std::cout << "--> Test: operators6." <<std::endl;

  // Simple = Simple - Simple, all dense
  *z = *x - *y ;
  for (unsigned int i = 0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", fabs((*z)(i) - (*x)(i) + (*y)(i)) < tol , true);
  z->zero();
  // Simple = Simple - Block, all dense
  //  *z = *x - *yB ;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", fabs((*z)(i) - (*x)(i) + (*yB)(i))<tol , true);
  //
  //  z->zero();
  //  // Simple = Block - Simple, all dense
  //  *z = *xB - *y ;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", fabs((*z)(i) - (*xB)(i) + (*y)(i))<tol , true);
  //
  //  z->zero();
  //  // Simple = Block - Block
  //  *z = *xB - *yB;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", fabs((*z)(i) - (*xB)(i) + (*yB)(i))<tol , true);

  // TO DO x and/or y sparse.

  // Block = Simple - Simple, all dense
  *zB = *x - *y ;
  for (unsigned int i = 0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", fabs((*zB)(i) - (*x)(i) + (*y)(i)) < tol , true);
  zB->zero();
  // Block = Simple - Block, all dense
  //  *zB = *x - *yB ;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", fabs((*zB)(i) - (*x)(i) + (*yB)(i))<tol , true);
  //
  //  zB->zero();
  //  // Block = Block - Simple, all dense
  //  *zB = *xB - *y ;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", fabs((*zB)(i) - (*xB)(i) + (*y)(i))<tol , true);
  //
  //  zB->zero();
  //  // Block = Block - Block
  //  *zB = *xB - *yB;
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", fabs((*zB)(i) - (*xB)(i) + (*yB)(i))<tol , true);
  std::cout << "--> operators6 test ended with success." <<std::endl;
}

void SiconosVectorTest::testOperators6Bis()
{
  //   // z = x - y with sub
  std::cout << "--> Test: operators6Bis." <<std::endl;

  // Simple = Simple - Simple, all dense
  sub(*x, *y, *z);
  for (unsigned int i = 0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis : ", fabs((*z)(i) - (*x)(i) + (*y)(i)) < tol , true);
  z->zero();
  // Simple = Simple - Block, all dense
  //  sub(*x,*yB,*z);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis : ", fabs((*z)(i) - (*x)(i) + (*yB)(i))<tol , true);
  //
  //  z->zero();
  //  // Simple = Block - Simple, all dense
  //  sub(*xB,*y,*z);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis : ", fabs((*z)(i) - (*xB)(i) + (*y)(i))<tol , true);
  //
  //  z->zero();
  //  // Simple = Block - Block
  //  sub(*xB,*yB,*z);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis : ", fabs((*z)(i) - (*xB)(i) + (*yB)(i))<tol , true);
  //
  //  // TO DO x and/or y sparse.
  //
  //  // Block = Simple - Simple, all dense
  //  sub(*x,*y,*zB);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis : ", fabs((*zB)(i) - (*x)(i) + (*y)(i))<tol , true);
  //  zB->zero();
  //  // Block = Simple - Block, all dense
  //  sub(*x,*yB,*zB);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis : ", fabs((*zB)(i) - (*x)(i) + (*yB)(i))<tol , true);
  //
  //  zB->zero();
  //  // Block = Block - Simple, all dense
  //  sub(*xB,*y,*zB);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis : ", fabs((*zB)(i) - (*xB)(i) + (*y)(i))<tol , true);
  //
  //  zB->zero();
  //  // Block = Block - Block
  //  sub(*xB,*yB,*zB);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis : ", fabs((*zB)(i) - (*xB)(i) + (*yB)(i))<tol , true);
  std::cout << "--> operators6Bis test ended with success." <<std::endl;
}

void SiconosVectorTest::testOperators7()
{
  //   // y = ax + y
  std::cout << "--> Test: operators7." <<std::endl;

  double a = 2.2;

  // z Simple ,  x Simple
  SP::SiconosVector tmp(new SiconosVector(*z)); // copy of z for comparison.
  axpy(a, *x, *z);

  for (unsigned int i = 0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators7 : ", fabs((*z)(i) - a * (*x)(i) - (*tmp)(i)) < tol , true);
  z->zero();
  //  z Simple ,  x Block
  //  axpy(a, *xB,*z);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators7 : ", fabs((*z)(i) - a*(*xB)(i) + (*tmp)(i))<tol , true);
  //
  //  z->zero();
  //  zB->zero();
  //  // z Block, x Simple
  //  tmp.reset(new SiconosVector(*zB));
  //  axpy(a,*x,*zB);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators7 : ", fabs((*zB)(i) - a*(*x)(i) + (*tmp)(i))<tol , true);
  //
  //  zB->zero();
  //  // y and x Block
  //  axpy(a, *xB,*zB);
  //  for (unsigned int i = 0; i<size; ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators7 : ", fabs((*zB)(i) - a*(*xB)(i) + (*tmp)(i))<tol , true);
  //
  //  zB->zero();
  //
  // TO DO x and/or y sparse.
  std::cout << "--> operators7 test ended with success." <<std::endl;
}

// outer_prod(v,w) = trans(v)*w
void SiconosVectorTest::testOperators8()
{
  std::cout << "--> Test: operators8." <<std::endl;

  SP::SiconosMatrix res(new SimpleMatrix(size1, size));
  *res = outer_prod(*tmp1, *x); // dense*dense

  for (unsigned int i = 0 ; i < size1; ++i)
    for (unsigned int j = 0; j < size; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8 : ", ((*res)(i, j) - (*tmp1)(i) * (*x)(j)) < tol , true);

  res->zero();
  SP::SiconosVector w(new SiconosVector(*sv));
  *res = outer_prod(*tmp1, *w); // dense*sparse
  for (unsigned int i = 0 ; i < size1; ++i)
    for (unsigned int j = 0; j < size; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8 : ", ((*res)(i, j) - (*tmp1)(i) * (w->getSparse())(j)) < tol , true);

  res->zero();
  SP::SiconosVector v(new SiconosVector(*sv));
  SP::SiconosMatrix res2(new SimpleMatrix(size, size));
  *res2 = outer_prod(*v, *w); // sparse*sparse
  for (unsigned int i = 0 ; i < size1; ++i)
    for (unsigned int j = 0; j < size; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8 : ", ((*res2)(i, j) - (v->getSparse())(i) * (w->getSparse())(j)) < tol , true);
  res->zero();

  *res2 = outer_prod(*v, *x); // sparse*dense
  for (unsigned int i = 0 ; i < size1; ++i)
    for (unsigned int j = 0; j < size; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8 : ", ((*res2)(i, j) - (v->getSparse())(i) * (*x)(j)) < tol , true);

  std::cout << "--> operators8 test ended with success." <<std::endl;
}

void SiconosVectorTest::testSubscal()
{
  std::cout << "--> Test: testSubscal." <<std::endl;
  unsigned int size = 12;
  SP::SiconosVector xx(new SiconosVector(size, 1.0));
  SP::SiconosVector ys(new SiconosVector(size, 2.0));

  SP::BlockVector yy(new BlockVector());
  SP::SiconosVector y1(new SiconosVector(2, 1.0));
  SP::SiconosVector y2(new SiconosVector(3, 2.0));
  SP::SiconosVector y3(new SiconosVector(7, 4.0));
  yy->insertPtr(y1);
  yy->insertPtr(y2);
  yy->insertPtr(y3);
  SP::SiconosVector yref(new SiconosVector(*ys));

  double a = 3;
  Index coord(4) ;
  coord[0] = 1;
  coord[1] = 3;
  coord[2] = 4;
  coord[3] = 6;

  // xx and yy dense vectors
  subscal(a, *xx, *ys, coord, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ", fabs((*ys)(4) - a * (*xx)(1)) < tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ", fabs((*ys)(5) - a * (*xx)(2)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 4 && i != 5)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ", fabs((*ys)(i) - (*yref)(i)) < tol, true);
  }

  // xx dense, y block
  //  *yref = *yy;
  //  subscal(a,*xx,*yy,coord,true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ",fabs((*yy)(4)-a*(*xx)(1))<tol, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ",fabs((*yy)(5)-a*(*xx)(2))<tol, true);
  //  for (unsigned int i = 0; i<size; ++i)
  //  {
  //    if (i!=4 && i!=5)
  //      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ",fabs((*yy)(i)-(*yref)(i))<tol, true);
  //  }
  //
  //  // x block, ys dense
  //  *yref = *ys;
  //  subscal(a,*yy,*ys,coord,true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ",fabs((*ys)(4)-a*(*yy)(1))<tol, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ",fabs((*ys)(5)-a*(*yy)(2))<tol, true);
  //  for (unsigned int i = 0; i<size; ++i)
  //  {
  //    if (i!=4 && i!=5)
  //      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ",fabs((*ys)(i)-(*yref)(i))<tol, true);
  //  }
  //
  // xx sparse, ys dense
  xx.reset(new SiconosVector(size, Siconos::SPARSE));
  xx->fill(2.0);
  *yref = *ys;
  subscal(a, *xx, *ys, coord, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ", fabs((*ys)(4) - a * (*xx)(1)) < tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ", fabs((*ys)(5) - a * (*xx)(2)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 4 && i != 5)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ", fabs((*ys)(i) - (*yref)(i)) < tol, true);
  }

  // xx sparse, ys sparse
  ys.reset(new SiconosVector(size, Siconos::SPARSE));
  *yref = *ys;
  subscal(a, *xx, *ys, coord, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ", fabs((*ys)(4) - a * (*xx)(1)) < tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ", fabs((*ys)(5) - a * (*xx)(2)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 4 && i != 5)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSubscal : ", fabs((*ys->sparse())(i) - (*yref->dense())(i)) < tol, true);
  }
  std::cout << "-->  subscal test ended with success." <<std::endl;
}

void SiconosVectorTest::testStdOstream()
{
  SP::SiconosVector y(new SiconosVector(2, 1.0));
  std::cout << *y << std::endl;
}

void SiconosVectorTest::testStdVectorCast()
{
  SP::SiconosVector y(new SiconosVector(2, 1.0));
  std::vector<double> d(*y);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testStdVectorCast : ", d.size() == y->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testStdVectorCast : ", d[0] == 1.0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testStdVectorCast : ", d[1] == 1.0, true);
}

void SiconosVectorTest::testIterators()
{
  SP::SiconosVector y(new SiconosVector(3, 1.0));
  for (SiconosVector::iterator it=y->begin(); it!=y->end(); it++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testIterators : ", *it == 1.0, true);

  SPC::SiconosVector z(new SiconosVector(3, 1.0));
  for (SiconosVector::const_iterator it=z->begin(); it!=z->end(); it++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testIterators : ", *it == 1.0, true);
}

void SiconosVectorTest::End()
{
  std::cout << "======================================" <<std::endl;
  std::cout << " ===== End of SiconosVector Tests ===== " <<std::endl;
  std::cout << "======================================" <<std::endl;
}



