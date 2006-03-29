/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
  nsv = new BlockVector(*u3);
  (static_cast<BlockVector*>(nsv))->add(*u4);

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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector : ", v->size() == 0, true);
  delete v;
  cout << " Constructor SV 0 ok" << endl;
}
*/
// from a file, see testRead

// Copy from a std vector
void SimpleVectorTest::testBuildSimpleVector1()
{
  cout << "====================================" << endl;
  cout << "=== Simple Vector tests start ...=== " << endl;
  cout << "====================================" << endl;
  cout << "--> Test: constructor 1." << endl;
  SimpleVector * v = new SimpleVector(vq);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector1 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector1 : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector1 : ", (v->getValues())(i) == vq[i], true);
  delete v;
  cout << "--> Constructor 1 test ended with success." << endl;
}

// Copy
void SimpleVectorTest::testBuildSimpleVector2()
{
  cout << "--> Test: constructor 2." << endl;
  SimpleVector *v = new SimpleVector(*u);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->size() == u->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", *v == *u, true);
  delete v;
  cout << "--> Constructor 2 test ended with success." << endl;
}

// copy from SiconosVector
void SimpleVectorTest::testBuildSimpleVector3()
{
  cout << "--> Test: constructor 3." << endl;
  // from a siconos which is a simple
  SiconosVector * tmp = new SimpleVector(vq);
  SimpleVector *v = new SimpleVector(*tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->size() == tmp->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", *v == *tmp, true);
  delete v;
  delete tmp;

  // from a siconos which is a block vector
  tmp = static_cast<SiconosVector*>(u2);

  SiconosVector * tmp2 = new BlockVector(*tmp);
  (static_cast<BlockVector*>(tmp2))->add(*u2);
  v = new SimpleVector(*tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", v->size() == 10, true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector2 : ", (*v)(i) == 1, true);
  delete v;
  delete tmp2;
  cout << "--> Constructor 3 test ended with success." << endl;
}

// with size
void SimpleVectorTest::testBuildSimpleVector4()
{
  cout << "--> Test: constructor 4." << endl;
  unsigned int SIZE = 10;
  SimpleVector *v = new SimpleVector(SIZE);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->size() == SIZE, true);
  delete v;
  cout << "--> Constructor 4 test ended with success." << endl;
}

// destructor

// zero
void SimpleVectorTest::testZero()
{
  cout << "--> Test: zero." << endl;
  SimpleVector *v = new SimpleVector(vq);
  v->zero();
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*v)(i) == 0, true);
  delete v;
  cout << "--> zero test ended with success." << endl;
}

// toString
// display

// setValue
void SimpleVectorTest::testSetValue()
{
  cout << "--> Test: setValue." << endl;
  double a = 4;
  int i = 2;
  u->setValue(i, a);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValue : ", (*u)(2) == 4, true);
  cout << "-->  setValue test ended with success." << endl;
}

// getValue
void SimpleVectorTest::testGetValue()
{
  cout << "--> Test: getValue." << endl;
  int i = 2;
  u->setValue(i, 8);
  u->getValue(i);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetValue : ", (*u)(2) == 8, true);
  cout << "-->  getValue test ended with success." << endl;
}

// setValues
void SimpleVectorTest::testSetValues()
{
  cout << "--> Test: setValues." << endl;
  u->setValues(vq);
  for (unsigned int i = 0; i < u->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValues : ", (u->getValues())(i) == vq[i], true);
  cout << "-->  setValues test ended with success." << endl;
}

// getValues
void SimpleVectorTest::testGetValues()
{
  cout << "--> Test: getValues." << endl;
  for (unsigned int i = 0; i < u->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetValues : ", (u2->getValues())(i) == 1, true);
  cout << "-->  getValues test ended with success." << endl;
}

// size
void SimpleVectorTest::testSize()
{
  cout << "--> Test: setSize." << endl;
  int i = u2->size();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSize : ", i == 5, true);
  cout << "-->  setSize test ended with success." << endl;
}

// write, read and read from a file constructor
void SimpleVectorTest::testReadWrite()
{
  cout << "--> Test: readWrite." << endl;
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
  cout << "--> readWrite test ended with success." << endl;
}

// getArray and norm: embedded blas or lapack functions => no tests

// norm

// OPERATORS

// += -=
void SimpleVectorTest::testOperatorPlusEqual()
{
  cout << "--> Test: operatorPlusEqual." << endl;
  SimpleVector *v = new SimpleVector(vq);

  SiconosVector *sv = new SimpleVector(vq);
  *v += *sv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == 2 * (*sv)(i), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isBlock(), false);

  *v -= *sv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == vq[i], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isBlock(), false);

  delete sv;

  *v += *nsv;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(2) == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(3) == 9, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(4) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isBlock(), false);

  v->zero();
  *v -= *nsv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", v->size() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(2) == -12, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(3) == -8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(4) == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isBlock(), false);
  delete v;
  cout << "--> operatorPlusEqual test ended with success." << endl;
}

// =
void SimpleVectorTest::testOperatorEqual()
{
  cout << "--> Test: operatorEqual." << endl;
  SimpleVector *v = new SimpleVector(vq);
  SimpleVector *w = new SimpleVector(vdotq);
  SiconosVector *z = new SimpleVector(vq);

  *v = *w ;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", v->size() == w->size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == (*w)(i), true);

  *v = *z;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", v->size() == z->size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == (*z)(i), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isBlock(), false);

  *v = *nsv;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorEqual : ", v->size() == nsv->size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorPlusEqualGEN : ", (*v)(i) == (*nsv)(i), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector3 : ", v->isBlock(), false);


  delete z;
  delete w;
  delete v;
  cout << "--> operatorEqual test ended with success." << endl;
}

// *= , /=

void SimpleVectorTest::testOperatorMultDivEqual()
{
  cout << "--> Test: operatorMultDivEqual." << endl;

  double d = 2;

  *u2 *= d;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualSPC : ", u2->size() == 5, true);
  for (unsigned int i = 0; i < u2->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorMultEqualGEN : ", (*u2)(i) == 2, true);

  *u2 /= d;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualSPC : ", u2->size() == 5, true);
  for (unsigned int i = 0; i < u2->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorDivEqualGEN : ", (*u2)(i) == 1, true);
  cout << "--> operatorMultDivEqual test ended with success." << endl;
}

// addition
void SimpleVectorTest::testAddition()
{
  cout << "--> Test: addition." << endl;
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
  cout << "--> addition test ended with success." << endl;
}

// subtraction
void SimpleVectorTest::testSubtraction()
{
  cout << "--> Test: subtraction ." << endl;
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
  cout << "--> subtraction test ended with success." << endl;
}
// +

void SimpleVectorTest::testExternalOperatorPlusMoins()
{
  cout << "--> Test: externalOperatorPlusMoins ." << endl;
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
  cout << "--> externalOperatorPlusMoins test ended with success." << endl;
}

// * /
void SimpleVectorTest::testExternalOperatorMultDiv()
{
  cout << "--> Test: externalOperatorMultDiv ." << endl;
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
  cout << "--> externalOperatorMultDiv test ended with success." << endl;
}

// matTransVectMult

void SimpleVectorTest::testExternalOperatorMultMat()
{
  cout << "--> Test: externalOperatorMultMat ." << endl;
  SimpleMatrix m(2, 4);
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

  SiconosVector * sv3 = new BlockVector(*u4);
  (static_cast<BlockVector*>(sv3))->add(*u4);
  (*sv3)(0) = 1;
  (*sv3)(1) = 2;
  (*sv3)(2) = 3;
  (*sv3)(3) = 4;
  sv = m * (*sv3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);


  delete sv3;
  delete sv2;
  delete v;
  cout << "--> externalOperatorMultMat test ended with success." << endl;
}

void SimpleVectorTest::testExternalOperatorMatTransMult()
{
  cout << "--> Test: externalOperatorMatTransMult ." << endl;
  SimpleMatrix m(4, 2);
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

  SiconosVector * sv3 = new BlockVector(*u4);
  (static_cast<BlockVector*>(sv3))->add(*u4);
  (*sv3)(0) = 1;
  (*sv3)(1) = 2;
  (*sv3)(2) = 3;
  (*sv3)(3) = 4;
  sv = matTransVecMult(m, *sv3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExternalOperatorMultMat : ", sv == res, true);


  delete sv3;
  delete sv2;
  delete v;
  cout << "--> externalOperatorMatTransMult test ended with success." << endl;
}

// ()
void SimpleVectorTest::testOperatorAccess()
{
  cout << "--> Test: operatorAccess ." << endl;
  SimpleVector *v = new SimpleVector(vq);

  (*v)(2) = 10;

  const double d = (*v)(1);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", v->size() == vq.size() , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", (*v)(2) != vq[2] , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", (v->getValues())(2) == 10 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperatorAccessRef : ", d == 1, true);
  delete v;
  cout << "--> operatorAccess test ended with success." << endl;
}

void SimpleVectorTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of SimpleVector Tests ===== " << endl;
  cout << "======================================" << endl;
}



