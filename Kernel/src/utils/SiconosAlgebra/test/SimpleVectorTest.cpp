/* Siconos version 1.0, Copyright INRIA 2005.
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

#include "SimpleVectorTest.h"
using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION(SimpleVectorTest);


void SimpleVectorTest::setUp()
{
  u = new SimpleVector(5);
  u->zero();
  vq.resize(5, 1);
  vdotq.resize(5, 1);
  for (int i = 0; i < 5; i++)
  {
    vq.at(i) = 1;
    vdotq.at(i) = 2;
  }

  u2 = new SimpleVector(vq);
  u3 = new SimpleVector(3);
  u4 = new SimpleVector(2);
  u3->zero();
  u4->zero();
  u3 ->setValue(2, 12);
  u4 ->setValue(0, 8);
  nsv = new CompositeVector(*u3);
  (static_cast<CompositeVector*>(nsv))->add(*u4);

}

void SimpleVectorTest::tearDown()
{
  delete nsv;
  delete u4;
  delete u3;
  delete u2;
  delete u;
}

// CONSTRUCTORS TESTS

// Default  ->private
/*void SimpleVectorTest::testBuildSimpleVector()
{
  SimpleVector * v = new SimpleVector();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector : ", v->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector : ", v->size() == 0, true);
  delete v;
  cout << " Constructor SV 0 ok" << endl;
}
*/
// from a file, see testRead

// Copy from a std vector
void SimpleVectorTest::testBuildSimpleVector1()
{
  SimpleVector * v = new SimpleVector(vq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector1 : ", v->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector1 : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector1 : ", (v->getValues())(i) == vq[i], true);
  cout << " Constructor SV 1 ok" << endl;
  delete v;
}

// Copy
void SimpleVectorTest::testBuildSimpleVector2()
{
  SimpleVector *v = new SimpleVector(*u);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->size() == u->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->getValues() == u->getValues(), true);
  cout << "SimpleVectorTest >>> testBuildSimpleVector2 ............................... OK\n ";
  delete v;
}

// copy from SiconosVector
void SimpleVectorTest::testBuildSimpleVector3()
{
  // from a siconos which is a simple
  SiconosVector * tmp = new SimpleVector(vq);
  SimpleVector *v = new SimpleVector(*tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->size() == tmp->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->getValues() == tmp->getValues(), true);
  cout << "SimpleVectorTest >>> testBuildSimpleVector2 ............................... OK\n ";
  delete v;
  delete tmp;

  // from a siconos which is a composite
  tmp = static_cast<SiconosVector*>(u2);

  SiconosVector * tmp2 = new CompositeVector(*tmp);
  (static_cast<CompositeVector*>(tmp2))->add(*u2);
  v = new SimpleVector(*tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->size() == 10, true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", (*v)(i) == 1, true);
  cout << "SimpleVectorTest >>> testBuildSimpleVector2 ............................... OK\n ";
  delete v;
  delete tmp2;

}

// with size
void SimpleVectorTest::testBuildSimpleVector4()
{
  unsigned int SIZE = 10;
  SimpleVector *v = new SimpleVector(SIZE);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->size() == SIZE, true);
  cout << "SimpleVectorTest >>> testBuildSimpleVector3 ............................... OK\n ";
  delete v;
}

// destructor

// zero
void SimpleVectorTest::testZero()
{
  SimpleVector *v = new SimpleVector(vq);
  v->zero();
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*v)(i) == 0, true);
  delete v;
}

// toString
// display

// setValue
void SimpleVectorTest::testSetValue()
{
  double a = 4;
  int i = 2;
  u->setValue(i, a);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValue : ", (*u)(2) == 4, true);
}

// getValue
void SimpleVectorTest::testGetValue()
{
  int i = 2;
  u->setValue(i, 8);
  u->getValue(i);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetValue : ", (*u)(2) == 8, true);
}

// setValues
void SimpleVectorTest::testSetValues()
{
  u->setValues(vq);
  for (unsigned int i = 0; i < u->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", (u->getValues())(i) == vq[i], true);
}

// getValues
void SimpleVectorTest::testGetValues()
{
  for (unsigned int i = 0; i < u->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetValues : ", (u2->getValues())(i) == 1, true);
}

// size
void SimpleVectorTest::testSize()
{
  int i = u2->size();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSize : ", i == 5, true);
}

// write, read and read from a file constructor
void SimpleVectorTest::testReadWrite()
{
  // write u2 into an ascii file
  bool isok = u2->write("testWrite_ascii.dat", "ascii");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testReadWrite : ", isok == true, true);

  // write u2 into a binary file
  bool isok2 = u2->write("testWrite_bin.dat", "binary");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testReadWrite : ", isok2 == true, true);

  // load v from an ascii file
  string input = "testWrite_ascii.dat";
  SimpleVector * v = new SimpleVector(input, true);
  // load v2 from a binary file
  string input_bin = "testWrite_bin.dat";
  SimpleVector *v2 = new SimpleVector(input_bin, false);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", v->size() == v2->size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", (*v)(i) == (*v2)(i), true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testRead : ", (*v)(i) == (*u2)(i), true);
  }
  delete v2;
  delete v;
  cout << "SimpleVectorTest >>> testWrite ............................... OK\n ";
}

// getArray and norm: embedded blas or lapack functions => no tests

// norm

// OPERATORS

// += -=
void SimpleVectorTest::testOperatorPlusEqual()
{
  SimpleVector *v = new SimpleVector(vq);

  SiconosVector *sv = new SimpleVector(vq);
  *v += *sv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == 2 * (*sv)(i), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isComposite(), false);

  *v -= *sv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == vq[i], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isComposite(), false);

  delete sv;

  *v += *nsv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(2) == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(3) == 9, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(4) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isComposite(), false);

  v->zero();
  *v -= *nsv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(2) == -12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(3) == -8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(4) == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isComposite(), false);
  delete v;
  cout << "SimpleVectorTest >>> testOperatorPlusEqualGEN ............................... OK\n ";
}

// =
void SimpleVectorTest::testOperatorEqual()
{
  SimpleVector *v = new SimpleVector(vq);
  SimpleVector *w = new SimpleVector(vdotq);
  SiconosVector *z = new SimpleVector(vq);

  *v = *w ;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", v->size() == w->size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == (*w)(i), true);

  *v = *z;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", v->size() == z->size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == (*z)(i), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isComposite(), false);

  *v = *nsv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", v->size() == nsv->size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == (*nsv)(i), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isComposite(), false);


  delete z;
  delete w;
  delete v;
  cout << "SimpleVectorTest >>> testOperatorEqual ............................... OK\n ";
}

// ==, !=

void SimpleVectorTest::testOperatorComp()
{
  SiconosVector *v = new SimpleVector(vq);
  SiconosVector *w = new SimpleVector(3);
  SiconosVector *z = new SimpleVector(vdotq);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", *v == *u2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", !(*v == *w), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", !(*v == *z), true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", !(*u == *u2), true);
  *u = *u2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", *u == *u2, true);


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", !(*v != *u2), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", *v != *w, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", *v != *z, true);

  delete z;
  delete w;
  delete v;
  cout << "SimpleVectorTest >>> testOperatorComp ............................... OK\n ";
}


// *= , /=

void SimpleVectorTest::testOperatorMultDivEqual()
{

  double d = 2;

  *u2 *= d;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualSPC : ", u2->size() == 5, true);
  for (unsigned int i = 0; i < u2->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*u2)(i) == 2, true);

  *u2 /= d;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", u2->size() == 5, true);
  for (unsigned int i = 0; i < u2->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualGEN : ", (*u2)(i) == 1, true);
  cout << "SimpleVectorTest >>> testOperatorMultDivEqual ............................... OK\n ";
}

// addition
void SimpleVectorTest::testAddition()
{
  SiconosVector * sv = new SimpleVector(*u2);

  *u = u2->addition(*sv);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", u->size() == 5, true);
  for (unsigned int i = 0; i < u->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*u)(i) == 2, true);

  *u = u2->addition(*nsv);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", u->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*u)(2) == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*u)(3) == 9, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*u)(4) == 1, true);
  delete sv;
  cout << "SimpleVectorTest >>> testAddition ............................... OK\n ";
}

// subtraction
void SimpleVectorTest::testSubtraction()
{
  SiconosVector * sv = new SimpleVector(*u2);

  *u = u2->subtraction(*sv);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", u->size() == 5, true);
  for (unsigned int i = 0; i < u->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*u)(i) == 0, true);

  *u = u2->subtraction(*nsv);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", u->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*u)(2) == -11, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*u)(3) == -7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*u)(4) == 1, true);

  delete sv;
  cout << "SimpleVectorTest >>> testSubtraction ............................... OK\n ";
}
// +

void SimpleVectorTest::testExternalOperatorPlusMoins()
{
  SimpleVector *v = new SimpleVector(5);
  SiconosVector *nsv1 = new SimpleVector(vq);

  *v = *nsv1 + *nsv1;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == 2, true);

  *v = *nsv1 - *nsv1;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == 0, true);

  delete nsv1;
  SimpleVector *v2 = new SimpleVector(vdotq);

  *v = *u2 + *v2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == 3, true);

  *v = *u2 - *v2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == -1, true);

  *v = *nsv + *nsv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(2) == 24, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(3) == 16, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(4) == 0, true);

  *v = *nsv - *nsv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(2) == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(3) == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(4) == 0, true);

  *v = *nsv + *v2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(2) == 14, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(3) == 10, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(4) == 2, true);

  *v = *nsv - *v2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(2) == 10, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(3) == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(4) == -2, true);

  delete v2;
  delete v;
  cout << "SimpleVectorTest >>> testExternalOperatorPlusMoins ............................... OK\n ";
}

// * /
void SimpleVectorTest::testExternalOperatorMultDiv()
{
  SimpleVector *v = new SimpleVector(5);
  double d = 2;

  *v = d * (*u2);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == 2, true);

  *v = (*u2) / d;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorPlusGEN : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == 0.5, true);

  delete v;
  cout << "SimpleVectorTest >>> testExternalOperatorMultDiv ............................... OK\n ";
}

// matTransVectMult

void SimpleVectorTest::testExternalOperatorMultMat()
{
  SiconosMatrix m(2, 4);
  m(0, 0) = 0;
  m(0, 1) = 1;
  m(0, 2) = -1;
  m(0, 3) = 0;
  m(1, 0) = 2;
  m(1, 1) = 1;
  m(1, 2) = -1;
  m(1, 3) = -2;

  SimpleVector *v = new SimpleVector(4);

  (*v)(0) = 1;
  (*v)(1) = 2;
  (*v)(2) = 3;
  (*v)(3) = 4;

  SimpleVector res(2);
  res(0) = -1;
  res(1) = -7;

  SimpleVector sv(2);
  sv = m * *v;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);

  SiconosVector * sv2 = new SimpleVector(4);
  (*sv2)(0) = 1;
  (*sv2)(1) = 2;
  (*sv2)(2) = 3;
  (*sv2)(3) = 4;
  sv = m * (*sv2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);

  SiconosVector * sv3 = new CompositeVector(*u4);
  (static_cast<CompositeVector*>(sv3))->add(*u4);
  (*sv3)(0) = 1;
  (*sv3)(1) = 2;
  (*sv3)(2) = 3;
  (*sv3)(3) = 4;
  sv = m * (*sv3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);


  delete sv3;
  delete sv2;
  delete v;
  cout << "SimpleVectorTest >>> testExternalOperatorMultMat ............................... OK\n ";
}

void SimpleVectorTest::testExternalOperatorMatTransMult()
{
  SiconosMatrix m(4, 2);
  m(0, 0) = 0;
  m(0, 1) = 2;
  m(1, 0) = 1;
  m(1, 1) = 1;
  m(2, 0) = -1;
  m(2, 1) = -1;
  m(3, 0) = 0;
  m(3, 1) = -2;

  SimpleVector *v = new SimpleVector(4);

  (*v)(0) = 1;
  (*v)(1) = 2;
  (*v)(2) = 3;
  (*v)(3) = 4;

  SimpleVector res(2);
  res(0) = -1;
  res(1) = -7;

  SimpleVector sv(2);
  sv = matTransVecMult(m, *v);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);

  SiconosVector * sv2 = new SimpleVector(4);
  (*sv2)(0) = 1;
  (*sv2)(1) = 2;
  (*sv2)(2) = 3;
  (*sv2)(3) = 4;
  sv = matTransVecMult(m, *sv2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);

  SiconosVector * sv3 = new CompositeVector(*u4);
  (static_cast<CompositeVector*>(sv3))->add(*u4);
  (*sv3)(0) = 1;
  (*sv3)(1) = 2;
  (*sv3)(2) = 3;
  (*sv3)(3) = 4;
  sv = matTransVecMult(m, *sv3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);


  delete sv3;
  delete sv2;
  delete v;
  cout << "SimpleVectorTest >>> testExternalOperatorMatTransMult ............................... OK\n ";
}

// ()
void SimpleVectorTest::testOperatorAccess()
{
  SimpleVector *v = new SimpleVector(vq);

  (*v)(2) = 10;

  const double d = (*v)(1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", v->isComposite(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", v->size() == vq.size() , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", (*v)(2) != vq[2] , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", (v->getValues())(2) == 10 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", d == 1, true);
  delete v;
  cout << "SimpleVectorTest >>> testOperatorAccess ............................... OK\n ";

}




