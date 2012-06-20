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
#include "BlockMatrixTest.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(BlockMatrixTest);


void BlockMatrixTest::setUp()
{
  tol = 1e-12;

  B.reset(new SimpleMatrix(2, 2, 1));
  C.reset(new SimpleMatrix(2, 4, 2));
  D.reset(new SimpleMatrix(2, 1, 3));
  E.reset(new SimpleMatrix(3, 2, 4));
  F.reset(new SimpleMatrix(3, 4, 5));
  G.reset(new SimpleMatrix(3, 1, 6));

  m.resize(6);
  m[0] = B ;
  m[1] = C ;
  m[2] = D ;
  m[3] = E ;
  m[4] = F ;
  m[5] = G ;

  tRow.resize(2);
  tRow[0] = 2;
  tRow[1] = 5;
  tCol.resize(3);
  tCol[0] = 2;
  tCol[1] = 6;
  tCol[2] = 7;

  mapRef.reset(new BlocksMat(2, 1));
  (*mapRef)(0, 0) = B;
  (*mapRef)(1, 0) = E;

}

void BlockMatrixTest::tearDown()
{
  m.clear();
}

//______________________________________________________________________________

void BlockMatrixTest::testConstructor0() // constructor with a vector of SP::SiconosMatrix
{
  cout << "====================================" << endl;
  cout << "=== Block Matrix tests start ...=== " << endl;
  cout << "====================================" << endl;
  cout << "--> Test: constructor 0." << endl;
  SP::SiconosMatrix test(new BlockMatrix(m, 2, 3));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size(0) == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size(1) == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(0, 1) == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(0, 2) == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(1, 0) == E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(1, 1) == F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(1, 2) == G, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", *test->tabRow() == tRow, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", *test->tabCol() == tCol, true);
  cout << "--> Constructor 0 test ended with success." << endl;
}

void BlockMatrixTest::testConstructor1() // Copy constructor, from a BlockMatrix
{
  cout << "--> Test: constructor 1." << endl;
  SP::BlockMatrix ref(new BlockMatrix(m, 2, 3));
  SP::SiconosMatrix test(new BlockMatrix(*ref));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->size(0) == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->size(1) == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(0, 0)) == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(0, 1)) == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(0, 2)) == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(1, 0)) == *E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(1, 1)) == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(1, 2)) == *G, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *test->tabRow() == tRow, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *test->tabCol() == tCol, true);
  cout << "--> Constructor 1(copy) test ended with success." << endl;
}

void BlockMatrixTest::testConstructor2() // Copy constructor, from a SiconosMatrix(Block)
{
  cout << "--> Test: constructor 2." << endl;
  SP::SiconosMatrix ref(new BlockMatrix(m, 2, 3));
  SP::SiconosMatrix test(new BlockMatrix(*ref));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->size(0) == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->size(1) == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(0, 0)) == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(0, 1)) == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(0, 2)) == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(1, 0)) == *E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(1, 1)) == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(1, 2)) == *G, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *test->tabRow() == tRow, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *test->tabCol() == tCol, true);

  cout << "--> Constructor 2(copy) test ended with success." << endl;
}
void BlockMatrixTest::testConstructor3() // Copy constructor, from a SiconosMatrix(Simple)
{
  cout << "--> Test: constructor 3." << endl;
  SP::SiconosMatrix ref(new SimpleMatrix(5, 7, 2.3));
  SP::SiconosMatrix test(new BlockMatrix(*ref));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->size(0) == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->size(1) == 7, true);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 7; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", fabs((*test)(i, j) - 2.3) < tol, true);
  SPC::Index tab = test->tabRow();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab->size() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[0] == 5, true);

  tab = test->tabCol();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab->size() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[0] == 7, true);

  cout << "--> Constructor 3(copy) test ended with success." << endl;
}

void BlockMatrixTest::testConstructor4() // Constructor from 4 SP::SiconosMatrix
{
  cout << "--> Test: constructor 4." << endl;
  SP::SiconosMatrix test(new BlockMatrix(B, C, E, F));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->size(0) == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->size(1) == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->block(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->block(0, 1) == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->block(1, 0) == E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->block(1, 1) == F, true);
  SPC::Index tab = test->tabRow();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab->size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[0] == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[1] == 5, true);
  tab = test->tabCol();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab->size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[0] == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[1] == 6, true);
  cout << "--> Constructor 4 test ended with success." << endl;
}

// void BlockMatrixTest::testResize()
// {
//   cout << "--> Test: resize." << endl;
//   SP::SiconosMatrix test(new BlockMatrix(m,2,3);
//   test->resize(3,4);
//   (*test)(2,0) = B;
//   (*test)(2,1) = C;
//   (*test)(2,2) = D;
//   (*test)(2,3) = C;
//   (*test)(0,3) = C;
//   (*test)(1,3) = F;

//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->isBlock() == true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->size(0) == 7, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->size(1) == 11, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(0,0) == B, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(0,1) == C, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(0,2) == D, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(0,3) == C, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(1,0) == E, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(1,1) == F, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(1,2) == G, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(1,3) == F, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(2,0) == B, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(2,1) == C, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(2,2) == D, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(2,3) == C, true);

//   cout << "--> resize test ended with success." << endl;
// }

void BlockMatrixTest::testNormInf()
{
  cout << "--> Test: normInf." << endl;
  SP::SiconosMatrix test(new BlockMatrix(m, 2, 3));
  test->zero();
  double n = 12;
  (*test)(4, 3) = n;
  (*test)(2, 1) = n - 3;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf: ", test->normInf() == n , true);
  cout << "--> normInf test ended with success." << endl;
}

void BlockMatrixTest::testZero()
{
  cout << "--> Test: zero." << endl;
  SP::SiconosMatrix A(new SimpleMatrix(2, 2));
  A->eye();
  SP::SiconosMatrix H(new SimpleMatrix(2, 4));
  H->eye();
  SP::SiconosMatrix I(new SimpleMatrix(5, 2));
  I->eye();
  SP::SiconosMatrix J(new SimpleMatrix(5, 4));
  J->eye();

  std::vector<SP::SiconosMatrix> v(4);
  v[0] = A ;
  v[1] = H ;
  v[2] = I ;
  v[3] = J ;
  SP::SiconosMatrix test(new BlockMatrix(v, 2, 2));
  test->zero();
  unsigned int n1 = test->size(0);
  unsigned int n2 = test->size(1);
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*test)(i, j) == 0, true);
  for (unsigned int i = 0; i < 2; ++i)
  {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*A)(i, j) == 0, true);
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*H)(i, j) == 0, true);
  }
  for (unsigned int i = 0; i < 5; ++i)
  {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*I)(i, j) == 0, true);
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*J)(i, j) == 0, true);
  }

  cout << "--> zero test ended with success." << endl;
}

void BlockMatrixTest::testEye()
{
  cout << "--> Test: eye." << endl;
  SP::SiconosMatrix A(new SimpleMatrix(2, 2));
  SP::SiconosMatrix H(new SimpleMatrix(2, 4));
  SP::SiconosMatrix I(new SimpleMatrix(5, 2));
  SP::SiconosMatrix J(new SimpleMatrix(5, 4));

  std::vector<SP::SiconosMatrix> v(4);
  v[0] = A ;
  v[1] = H ;
  v[2] = I ;
  v[3] = J ;
  SP::SiconosMatrix test(new BlockMatrix(v, 2, 2));
  test->eye();
  unsigned int n1 = test->size(0);
  unsigned int n2 = test->size(1);
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*test)(i, j) == 1, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*test)(i, j) == 0, true);

  for (unsigned int i = 0; i < 2; ++i)
  {
    for (unsigned int j = 0; j < 2; ++j)
    {
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*A)(i, j) == 1, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*A)(i, j) == 0, true);
    }
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*H)(i, j) == 0, true);
  }
  for (unsigned int i = 0; i < 5; ++i)
  {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*I)(i, j) == 0, true);
    for (unsigned int j = 0; j < 4; ++j)
    {
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*J)(i, j) == 1, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*J)(i, j) == 0, true);
    }
  }

  cout << "--> eye test ended with success." << endl;
}
// Add tests with getDense ...

void BlockMatrixTest::testGetSetRowCol()
{
  cout << "--> Test: get, set Row and Col." << endl;
  SP::SiconosVector tmp(new SiconosVector(6));
  SP::SiconosVector tmp1(new SiconosVector(6));
  (*tmp1)(0) = 1;
  (*tmp1)(2) = 2;
  SP::SiconosMatrix A(new SimpleMatrix(2, 2));
  SP::SiconosMatrix H(new SimpleMatrix(2, 4));
  SP::SiconosMatrix I(new SimpleMatrix(5, 2));
  SP::SiconosMatrix J(new SimpleMatrix(5, 4));
  std::vector<SP::SiconosMatrix> v(4);
  v[0] = A ;
  v[1] = H ;
  v[2] = I ;
  v[3] = J ;
  SP::SiconosMatrix test(new BlockMatrix(v, 2, 2));

  test->setRow(1, *tmp1);
  test->getRow(1, *tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", *tmp == *tmp1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*A)(1, 0) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*H)(1, 0) == 2, true);

  SP::SiconosVector tmp2(new SiconosVector(7));
  SP::SiconosVector tmp3(new SiconosVector(7));
  (*tmp3)(0) = 1;
  (*tmp3)(2) = 2;
  test->setCol(1, *tmp3);
  test->getCol(1, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", *tmp2 == *tmp3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*A)(0, 1) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*I)(0, 1) == 2, true);

  cout << "--> get, set Row and Col tests ended with success." << endl;
}

void BlockMatrixTest::testAssignment()
{
  cout << "--> Test: assignment." << endl;
  SP::SiconosMatrix Btmp(new SimpleMatrix(2, 2));
  SP::SiconosMatrix Ctmp(new SimpleMatrix(2, 5));
  SP::SiconosMatrix Dtmp(new SimpleMatrix(3, 2));
  SP::SiconosMatrix Etmp(new SimpleMatrix(3, 5));

  SP::SiconosMatrix test(new BlockMatrix(Btmp, Ctmp, Dtmp, Etmp));
  // Block = Siconos(Simple)
  unsigned int size0 = test->size(0), size1 = test->size(1);
  SP::SiconosMatrix ref(new SimpleMatrix(size0, size1));
  for (unsigned int i = 0; i < size0 ; ++i)
    for (unsigned int j = 0; j < size1; ++j)
      (*ref)(i, j) = i + j;
  *test = *ref;
  for (unsigned int i = 0; i < size0 ; ++i)
    for (unsigned int j = 0; j < size1; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", fabs((*test)(i, j) - (*ref)(i, j)) < tol , true);

  // Block = Siconos(Block)
  test->zero();

  ref.reset(new BlockMatrix(m, 2, 3));
  *test = *ref;
  for (unsigned int i = 0; i < size0 ; ++i)
    for (unsigned int j = 0; j < size1; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", fabs((*test)(i, j) - (*ref)(i, j)) < tol , true);

  // Block = Block
  test->zero();
  SP::BlockMatrix ref2(new BlockMatrix(m, 2, 3));
  *test = *ref2;
  for (unsigned int i = 0; i < size0 ; ++i)
    for (unsigned int j = 0; j < size1; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", fabs((*test)(i, j) - (*ref2)(i, j)) < tol , true);
  cout << "-->  test assignment ended with success." << endl;
}

void BlockMatrixTest::testOperators1()
{
  cout << "--> Test: operators1." << endl;
  double tol = 1e-10;
  SP::SiconosMatrix Ab(new BlockMatrix(m, 2, 3));
  SP::SiconosMatrix Cb(new BlockMatrix(*Ab));
  SP::SiconosMatrix A(new SimpleMatrix(5, 7));

  for (unsigned int i = 0; i < 5 ; ++i)
    for (unsigned int j = 0; j < 7; ++j)
      (*A)(i, j) = i + j;

  double a = 2.3;
  int a1 = 2;

  // Block *= scal or /= scal
  *Cb *= a;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - a * (*Ab)(i, j)) < tol , true);

  *Cb *= a1;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - a1 * a * (*Ab)(i, j)) < tol , true);

  *Cb /= a;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - a1 * (*Ab)(i, j)) < tol , true);
  *Cb /= a1;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - (*Ab)(i, j)) < tol , true);

  // Block +=  Simple
  *Cb += *A;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*A)(i, j)) < tol , true);

  // Block -=  Block
  *Cb -= *Ab;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - (*A)(i, j)) < tol , true);

  // Block += Block
  *Cb += *Ab;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*A)(i, j)) < tol , true);

  // Block -= Simple
  *Cb -= *A;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - (*Ab)(i, j)) < tol , true);

  cout << "-->  test operators1 ended with success." << endl;
}

void BlockMatrixTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of BlockMatrix Tests ===== " << endl;
  cout << "======================================" << endl;
}
