/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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

#include "MySimpleVectorTest.h"

using namespace std;
using namespace boost::numeric::ublas;

CPPUNIT_TEST_SUITE_REGISTRATION(MySimpleVectorTest);


void MySimpleVectorTest::setUp()
{
  tol = 1e-14;
  ref = new MySimpleVector(5);
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

void MySimpleVectorTest::tearDown()
{
  delete ref;
  delete dv;
  delete sv;
}

// Copy from a std vector
void MySimpleVectorTest::testConstructor1()
{
  cout << "====================================" << endl;
  cout << "=== My Simple Vector tests start ...=== " << endl;
  cout << "====================================" << endl;
  cout << "--> Test: constructor 1." << endl;
  MySimpleVector * v = new MySimpleVector(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->getNum() == 1, true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", (*v)(i) == 0, true);
  delete v;
  v = new MySimpleVector(3, SPARSE);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->getNum() == 4, true);
  delete v;
  cout << "--> Constructor 1 test ended with success." << endl;
}

void MySimpleVectorTest::testConstructor2()
{
  cout << "--> Test: constructor 2." << endl;

  MySimpleVector *v = new MySimpleVector(vq);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", v->size() == vq.size(), true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", (*v)(i) == vq[i], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", v->getNum() == 1, true);
  delete v;
  cout << "--> Constructor 2 test ended with success." << endl;
}

// copy from a MySiconosVector
void MySimpleVectorTest::testConstructor3()
{
  cout << "--> Test: constructor 3." << endl;
  MySiconosVector * tmp = new MySimpleVector(vq);
  MySimpleVector *v = new MySimpleVector(*tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", *v == *tmp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->getNum() == 1, true);
  delete v;
  delete tmp;
  tmp = new MySimpleVector(3, SPARSE);
  v = new MySimpleVector(*tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->size() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", *v == *tmp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", v->getNum() == 4, true);
  delete v;
  delete tmp;
  cout << "--> Constructor 3 test ended with success." << endl;
}

// with size
void MySimpleVectorTest::testConstructor4()
{
  cout << "--> Test: constructor 4." << endl;
  MySimpleVector * tmp = new MySimpleVector(vq);
  MySimpleVector *v = new MySimpleVector(*tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->size() == vq.size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", *v == *tmp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->getNum() == 1, true);
  delete v;
  delete tmp;
  tmp = new MySimpleVector(4, SPARSE);
  v = new MySimpleVector(*tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->size() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", *v == *tmp, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", v->getNum() == 4, true);
  delete v;
  delete tmp;
  cout << "--> Constructor 4 test ended with success." << endl;
}

void MySimpleVectorTest::testConstructor5()
{
  cout << "--> Test: constructor 5." << endl;
  MySiconosVector * v = new MySimpleVector(*dv);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->size() == dv->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", norm_inf(v->getDense() - *dv) < tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", v->getNum() == 1, true);
  delete v;
  cout << "--> Constructor 5 test ended with success." << endl;
}
void MySimpleVectorTest::testConstructor6()
{
  cout << "--> Test: constructor 6." << endl;
  MySiconosVector * v = new MySimpleVector(*sv);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", v->size() == sv->size(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", norm_inf(v->getSparse() - *sv) < tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", v->getNum() == 4, true);
  delete v;
  cout << "--> Constructor 6 test ended with success." << endl;
}

void MySimpleVectorTest::testConstructor7()
{
  cout << "--> Test: constructor 7." << endl;
  MySiconosVector * v = new MySimpleVector("vect.dat", 1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", v->size() == 4, true);
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", (*v)(i) == i, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", v->getNum() == 1, true);
  delete v;
  cout << "--> Constructor 7 test ended with success." << endl;
}
// destructor

// zero
void MySimpleVectorTest::testZero()
{
  cout << "--> Test: zero." << endl;
  MySiconosVector *v = new MySimpleVector(vq);
  v->zero();
  for (unsigned int i = 0; i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*v)(i) == 0, true);
  delete v;
  cout << "--> zero test ended with success." << endl;
}

void MySimpleVectorTest::testNorm()
{
  cout << "--> Test: norm." << endl;
  MySiconosVector *v = new MySimpleVector(*dv);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNorm : ", v->normInf() == norm_inf(*dv), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNorm : ", v->norm() == norm_2(*dv), true);
  delete v;
  cout << "--> norm test ended with success." << endl;
}

//resize
void MySimpleVectorTest::testResize()
{
  cout << "--> Test: resize." << endl;
  MySiconosVector *v = new MySimpleVector(vq);
  v->resize(6);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test resize : ", v->isBlock(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test resize : ", v->size() == 6, true);
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test resize : ", (*v)(i) == vq[i], true);
  for (unsigned int i = vq.size(); i < v->size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test resize : ", (*v)(i) == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test resize : ", v->getNum() == 1, true);
  delete v;
  cout << "-->  resize test ended with success." << endl;
}

// OPERATORS

// =
void MySimpleVectorTest::testAssignment()
{
  cout << "--> Test: operators1." << endl;
  MySiconosVector *v = new MySimpleVector(vq.size());
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == 0, true);
  MySiconosVector *w = new MySimpleVector(vq);
  *v = *w;
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment : ", (*v)(i) == (*w)(i), true);

  cout << "--> operators1 test ended with success." << endl;
}

// ()
void MySimpleVectorTest::testOperators1()
{
  cout << "--> Test: operators1." << endl;

  cout << "--> operators1 test ended with success." << endl;
}

// +=, -=, *=, /=
void MySimpleVectorTest::testOperators2()
{
  cout << "--> Test: operators2." << endl;
  MySiconosVector *v = new MySimpleVector(vq.size());
  MySiconosVector *w = new MySimpleVector(vq);

  // dense +=, -=, ... dense
  *v += *w;
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", (*v)(i) == (*w)(i), true);
  double mult = 2.2;
  int mult1 = 2;
  *v *= mult;
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", ((*v)(i) - mult * (*w)(i)) < tol, true);
  *v *= mult1;
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", ((*v)(i) - mult1 * mult * (*w)(i)) < tol, true);
  *v /= mult1;
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", ((*v)(i) - mult * (*w)(i)) < tol, true);
  *v /= mult;
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", ((*v)(i) - (*w)(i)) < tol, true);
  *v *= 3;
  *v -= *w;
  for (unsigned int i = 0; i < vq.size(); i++)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", ((*v)(i) - 2 * (*w)(i)) < tol, true);
  delete w;

  // dense +=, -=, ... sparse
  w = new MySimpleVector(*sv);
  v->resize(5);
  v->zero();
  *v += *w;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", norm_inf(v->getDense() - *sv) < tol, true);
  *v *= 3;
  *v -= *w;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", norm_inf(v->getDense() - 2 * *sv) < tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", v->getNum() == 1, true);
  delete v;
  // sparse += -= sparse
  v = new MySimpleVector(5, SPARSE);
  *v += *w;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", norm_inf(v->getSparse() - *sv) < tol, true);
  *v *= 3;
  *v -= *w;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", norm_inf(v->getSparse() - 2 * *sv) < tol, true);
  delete w;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators2 : ", v->getNum() == 4, true);
  delete v;
  cout << "--> operators2 test ended with success." << endl;
}

// inner_prod
void MySimpleVectorTest::testOperators3()
{
  cout << "--> Test: operators3." << endl;
  MySiconosVector *v = new MySimpleVector(vq);
  MySiconosVector *w = new MySimpleVector(vq);
  double res;
  // *
  res = inner_prod(*v, *w); // dense*dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", res == 55, true);
  delete w;
  w = new MySimpleVector(*sv);
  res = inner_prod(*v, *w); // dense*sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", res == 44, true);
  delete v;
  v = new MySimpleVector(*sv);
  res = inner_prod(*v, *w); // sparse*sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", res == 484, true);
  delete w;
  w = new MySimpleVector(vq);
  res = inner_prod(*v, *w); // sparse*dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", res == 44, true);
  cout << "--> operators3 test ended with success." << endl;
}

// vector * or / by double or int
void MySimpleVectorTest::testOperators4()
{
  cout << "--> Test: operators4." << endl;
  MySiconosVector *v = new MySimpleVector(vq);
  MySiconosVector * res = new MySimpleVector(5);
  double mult = 2.2;
  int mult1 = 2;
  // *
  *res = mult**v; // dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getDense() - mult * v->getDense()) < tol, true);
  res->zero();
  *res = *v * mult; // dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getDense() - mult * v->getDense()) < tol, true);
  res->zero();
  *res = mult1**v; // dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getDense() - mult1 * v->getDense()) < tol, true);
  res->zero();
  *res = *v * mult1; // dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getDense() - mult1 * v->getDense()) < tol, true);
  res->zero();
  *res = *v / mult; // dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getDense() - v->getDense() / mult) < tol, true);
  res->zero();
  *res = *v / mult1; // dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getDense() - v->getDense() / mult1) < tol, true);
  res->zero();
  delete v;
  v = new MySimpleVector(*sv);
  delete res;
  res = new MySimpleVector(5, SPARSE);
  *res = mult**v; // sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getSparse() - mult * v->getSparse()) < tol, true);
  res->zero();
  *res = *v * mult; // sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getSparse() - mult * v->getSparse()) < tol, true);
  res->zero();
  *res = mult1**v; // sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getSparse() - mult1 * v->getSparse()) < tol, true);
  res->zero();
  *res = *v * mult1; // sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getSparse() - mult1 * v->getSparse()) < tol, true);
  res->zero();
  *res = *v / mult; // sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getSparse() - v->getSparse() / mult) < tol, true);
  res->zero();
  *res = *v / mult1; // sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators4 : ", norm_2(res->getSparse() - v->getSparse() / mult1) < tol, true);
  delete res;
  delete v;

  cout << "--> operators4 test ended with success." << endl;
}

// +
void MySimpleVectorTest::testOperators5()
{
  cout << "--> Test: operators5." << endl;
  MySiconosVector *v = new MySimpleVector(vq);
  MySiconosVector *w = new MySimpleVector(vq);
  MySiconosVector * res = new MySimpleVector(5);
  MySiconosVector * res2 = new MySimpleVector(5, SPARSE);

  *res = *v + *w; // dense+dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", norm_2(res->getDense() - 2 * v->getDense()) < tol, true);
  delete w;
  w = new MySimpleVector(*sv);
  *res = *v + *w; // dense+sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", norm_2(res->getDense() - v->getDense() - w->getSparse()) < tol, true);
  delete v;
  v = new MySimpleVector(*sv);
  *res2 = *v + *w; // sparse+sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", norm_inf(res2->getSparse() - 2 * *sv) < tol, true);
  delete w;
  w = new MySimpleVector(vq);
  *res = *v + *w; // sparse+dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators5 : ", norm_2(res->getDense() - v->getSparse() - w->getDense()) < tol, true);
  delete v;
  delete w;
  delete res;
  delete res2;
  cout << "--> operators5 test ended with success." << endl;
}

// -
void MySimpleVectorTest::testOperators6()
{
  cout << "--> Test: operators6." << endl;
  MySiconosVector *v = new MySimpleVector(vq);
  MySiconosVector *w = new MySimpleVector(vq);
  MySiconosVector * res = new MySimpleVector(5);
  MySiconosVector * res2 = new MySimpleVector(5, SPARSE);

  *res = 2 * *v - *w; // dense-dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", norm_2(res->getDense() - v->getDense()) < tol, true);
  delete w;
  w = new MySimpleVector(*sv);
  *res = *v - *w; // dense-sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", norm_2(res->getDense() - v->getDense() + w->getSparse()) < tol, true);
  delete v;
  v = new MySimpleVector(*sv);
  *res2 = 2 * *v - *w; // sparse-sparse
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", norm_inf(res2->getSparse() - *sv) < tol, true);
  delete w;
  w = new MySimpleVector(vq);
  *res = *v - *w; // sparse-dense
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6 : ", norm_2(res->getDense() - v->getSparse() + w->getDense()) < tol, true);
  delete v;
  delete w;
  delete res;
  delete res2;
  cout << "--> operators6 test ended with success." << endl;
}

// outer_prod(v,w) = trans(v)*w
void MySimpleVectorTest::testOperators7()
{
  cout << "--> Test: operators7." << endl;
  MySiconosVector *v = new MySimpleVector(vq);
  v->resize(3);
  MySiconosVector *w = new MySimpleVector(vq);
  MySiconosMatrix *res = new MySimpleMatrix(3, 5);
  *res = outer_prod(*v, *w); // dense*dense
  for (unsigned int i = 0 ; i < 3; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", ((*res)(i, j) - (*v)(i) * (*w)(j)) < tol , true);
  delete w;
  res->zero();
  w = new MySimpleVector(*sv);
  *res = outer_prod(*v, *w); // dense*sparse
  for (unsigned int i = 0 ; i < 3; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", ((*res)(i, j) - (*v)(i) * (w->getSparse())(j)) < tol , true);
  delete v;
  res->zero();
  v = new MySimpleVector(*sv);
  *res = outer_prod(*v, *w); // sparse*sparse
  for (unsigned int i = 0 ; i < 3; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", ((*res)(i, j) - (v->getSparse())(i) * (w->getSparse())(j)) < tol , true);
  delete w;
  res->zero();
  w = new MySimpleVector(vq);
  *res = outer_prod(*v, *w); // sparse*dense
  for (unsigned int i = 0 ; i < 3; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators3 : ", ((*res)(i, j) - (v->getSparse())(i) * (*w)(j)) < tol , true);

  delete w;
  delete v;
  delete res;
  cout << "--> operators7 test ended with success." << endl;
}


void MySimpleVectorTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of MySimpleVector Tests ===== " << endl;
  cout << "======================================" << endl;
}



