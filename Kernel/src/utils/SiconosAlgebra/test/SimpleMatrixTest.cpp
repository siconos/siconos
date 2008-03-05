/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "SimpleMatrixTest.h"
#include "SimpleVector.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

CPPUNIT_TEST_SUITE_REGISTRATION(SimpleMatrixTest);

// Note FP: tests are (rather) complete for Dense objects but many are missing for other cases (Triang, Symm etc ...).

void SimpleMatrixTest::setUp()
{
  tol = 1e-9;

  fic1 = "mat1.dat"; // 2 X 2
  fic2 = "mat2.dat"; // 2 X 3
  SicM = new SimpleMatrix(fic1, 1);
  SimM = new SimpleMatrix(fic2, 1);

  std::vector<double> v3(2, 0);
  std::vector<double> v4(2, 0);
  std::vector<double> v5(3, 0);
  v4[0] = 6;
  v4[1] = 9;
  v5[0] = 8;
  v5[1] = 9;
  v5[2] = 10;

  vect1 = new SimpleVector(v3);
  vect2 = new SimpleVector(v4); // vect2 != vect1, but vect2 == SimM second column
  vect3 = new SimpleVector(v5); // vect3 != vect1, but vect3 == SimM second row

  // Dense
  D = new DenseMat(2, 2);
  for (unsigned i = 0; i < D->size1(); ++ i)
    for (unsigned j = 0; j < D->size2(); ++ j)
      (*D)(i, j) = 3 * i + j;
  // Triang
  T = new TriangMat(3, 3);
  for (unsigned i = 0; i < T->size1(); ++ i)
    for (unsigned j = i; j < T->size2(); ++ j)
      (*T)(i, j) = 3 * i + j;
  T2 = new TriangMat(4, 4);
  for (unsigned i = 0; i < T2->size1(); ++ i)
    for (unsigned j = i; j < T2->size2(); ++ j)
      (*T2)(i, j) = 3 * i + j;
  // Sym
  S = new SymMat(3, 3);
  for (unsigned i = 0; i < S->size1(); ++ i)
    for (unsigned j = i; j < S->size2(); ++ j)
      (*S)(i, j) = 3 * i + j;
  S2 = new SymMat(4, 4);
  for (unsigned i = 0; i < S2->size1(); ++ i)
    for (unsigned j = i; j < S2->size2(); ++ j)
      (*S2)(i, j) = 3 * i + j;
  // Sparse
  SP = new SparseMat(4, 4);
  for (unsigned i = 0; i < SP->size1(); ++ i)
    for (unsigned j = 0; j < SP->size2(); ++ j)
      (*SP)(i, j) = 3 * i + j;

  // Banded
  Band = new BandedMat(4, 4, 1, 1);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      (*Band)(i, j) = 3 * i + j;
  Band2 = new BandedMat(4, 3, 1, 1);
  for (signed i = 0; i < signed(Band2->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band2->size2())); ++ j)
      (*Band2)(i, j) = 3 * i + j;

  // Zero
  Z = new ZeroMat(3, 3);
  Z2 = new ZeroMat(4, 4);
  // Identity
  I = new IdentityMat(3, 3);
  I2 = new IdentityMat(4, 4);

  // BlockMat
  size = 10;
  size2 = 10;

  C = new SimpleMatrix(size, size);
  A = new SimpleMatrix("A.dat");
  B = new SimpleMatrix("B.dat");

  m1 = new SimpleMatrix(size - 2, size - 2, 1);
  m2 = new SimpleMatrix(size - 2, 2, 2);
  m3 = new SimpleMatrix(2, size - 2, 3);
  m4 = new SimpleMatrix(2, 2, 4);
  m5 = new SimpleMatrix(size2 - 2, size2 - 2, 1);
  m6 = new SimpleMatrix(size2 - 2, 2, 2);
  m7 = new SimpleMatrix(2, size2 - 2, 3);
  m8 = new SimpleMatrix(2, 2, 4);
  Ab = new BlockMatrix(m1, m2, m3, m4);
  Bb = new BlockMatrix(3 * *Ab);
  Cb = new BlockMatrix(m5, m6, m7, m8);


}

void SimpleMatrixTest::tearDown()
{

  delete Ab;
  delete Bb;
  delete Cb;
  delete Cb2;
  delete m1;
  delete m2;
  delete m3;
  delete m4;
  delete A;
  delete B;
  delete C;
  delete I;
  delete Z;
  delete I2;
  delete Z2;
  delete vect1;
  delete vect2;
  delete vect3;
  delete SicM;
  delete SimM;
  delete T;
  delete T2;
  delete SP;
  delete Band;
  delete Band2;

  delete D;
  delete S;
  delete S2;

}

//______________________________________________________________________________

void SimpleMatrixTest::testConstructor0() // constructor with TYP and dim
{
  cout << "====================================" << endl;
  cout << "=== Simple Matrix tests start ...=== " << endl;
  cout << "====================================" << endl;
  cout << "--> Test: constructor 0." << endl;
  SimpleMatrix * test = new SimpleMatrix(2, 3);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getNum() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size(0) == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size(1) == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->normInf() < tol, true);
  delete test;
  cout << "--> Constructor 0 test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor1() // Copy constructor, from a SimpleMatrix
{
  cout << "--> Test: constructor 1." << endl;
  SimpleMatrix * test = new SimpleMatrix(*SimM);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *test == *SimM, true);
  delete test;
  cout << "--> Constructor 1 (copy) test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor2() // Copy constructor, from a SiconosMatrix
{
  cout << "--> Test: constructor 2." << endl;
  SimpleMatrix *  test = new SimpleMatrix(*SicM);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *test == *SicM, true);
  delete test;
  cout << "--> Constructor 2 (copy) test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor3() // Copy constructor, from a BlockMatrix
{
  cout << "--> Test: constructor 3." << endl;
  SimpleMatrix *  test = new SimpleMatrix(*Ab);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", *test == *Ab, true);
  delete test;
  cout << "--> Constructor 3 (copy) test ended with success." << endl;
}

// void SimpleMatrixTest::testConstructor3()
// {
//   cout << "--> Test: constructor 3." << endl;
//   SimpleMatrix * test = new SimpleMatrix(SPARSE);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ",test->getNum() == 4, true);
//   delete test;
//   cout << "--> Constructor 3 test ended with success." << endl;
// }

void SimpleMatrixTest::testConstructor4()
{
  cout << "--> Test: constructor 4." << endl;
  SiconosMatrix * test = new SimpleMatrix(*D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->getNum() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", norm_inf(test->getDense() - *D) == 0, true);
  delete test;
  cout << "--> Constructor 4 test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor5()
{
  cout << "--> Test: constructor 5." << endl;
  SiconosMatrix * test = new SimpleMatrix(*T);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", test->getNum() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", norm_inf(test->getTriang() - *T) == 0, true);
  delete test;
  cout << "--> Constructor 5 test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor6()
{
  cout << "--> Test: constructor 6." << endl;
  SiconosMatrix * test = new SimpleMatrix(*S);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", test->getNum() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", norm_inf(test->getSym() - *S) == 0, true);
  delete test;
  cout << "--> Constructor 6 test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor7()
{
  cout << "--> Test: constructor 7." << endl;
  SiconosMatrix * test = new SimpleMatrix(*SP);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", test->getNum() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", norm_inf(test->getSparse() - *SP) == 0, true);
  delete test;
  cout << "--> Constructor 7 test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor8()
{
  cout << "--> Test: constructor 8." << endl;
  cout << "--> Constructor 8 test ended with success." << endl;
  SiconosMatrix * test = new SimpleMatrix(*Band);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor8 : ", test->getNum() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor8 : ", norm_inf(test->getBanded() - *Band) == 0, true);
  delete test;
}

void SimpleMatrixTest::testConstructor9() // constructor with TYP and dim and input value
{
  cout << "--> Test: constructor 9." << endl;
  SimpleMatrix * test = new SimpleMatrix(2, 3, 4.5);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test->getNum() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test->size(0) == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test->size(1) == 3, true);
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0 ; j < 3; ++j)
    {
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", (*test)(i, j) == 4.5, true);
    }
  delete test;
  cout << "--> Constructor 9 test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor10()
{
  cout << "--> Test: constructor 10." << endl;
  SiconosMatrix * test = new SimpleMatrix(fic1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor10 : ", *test == *SicM, true);
  delete test;
  cout << "--> Constructor 10 test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor11()
{
  cout << "--> Test: constructor 11." << endl;
  cout << "--> Constructor 11 test ended with success." << endl;
  SiconosMatrix * test = new SimpleMatrix(*Z);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor11 : ", test->getNum() == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor11 : ", test->normInf() == 0, true);
  delete test;
}

void SimpleMatrixTest::testConstructor12()
{
  cout << "--> Test: constructor 12." << endl;
  cout << "--> Constructor 12 test ended with success." << endl;
  SiconosMatrix * test = new SimpleMatrix(*I);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor12 : ", test->getNum() == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor12 : ", test->normInf() == 1, true);
  delete test;
}

// Add tests with getDense ...

void SimpleMatrixTest::testZero()
{
  cout << "--> Test: zero." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*SimM);
  tmp->zero();
  unsigned int n1 = tmp->size(0);
  unsigned int n2 = tmp->size(1);
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*tmp)(i, j) == 0, true);
  delete tmp;
  cout << "--> zero test ended with success." << endl;
}

void SimpleMatrixTest::testEye()
{
  cout << "--> Test: eye." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*SimM);
  tmp->eye();
  unsigned int n1 = tmp->size(0);
  unsigned int n2 = tmp->size(1);
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      if (i != j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*tmp)(i, j) == 0, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*tmp)(i, j) == 1, true);
  delete tmp;
  cout << "--> eye test ended with success." << endl;
}

void SimpleMatrixTest::testResize()
{
  cout << "--> Test: resize." << endl;
  SiconosMatrix * tmp = new SimpleMatrix(*SicM);
  tmp->resize(3, 4);
  unsigned int n1 = SicM->size(0);
  unsigned int n2 = SicM->size(1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(0) == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(1) == 4, true);

  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", fabs((*tmp)(i, j) - (*SicM)(i, j)) < tol, true);
  //   for(unsigned int i = n1; i<3; ++i)
  //     for(unsigned int j=0;j<4;++j)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", fabs((*tmp)(i,j)) < tol, true);
  //   for(unsigned int j = n2; j<4; ++j)
  //     for(unsigned int i=0;i<3;++i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", fabs((*tmp)(i,j)) < tol, true)
  ;
  // Check the effect of bool = false (ie preserve == false in boost resize)
  //   tmp->resize(6,8, false);
  //   tmp->display();
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(0) == 6, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(1) == 8, true);
  //   for(unsigned int i = 0; i<6; ++i)
  //     for(unsigned int j=0;j<8;++j)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", (*tmp)(i,j) == 0 , true);
  //   // Reduction ...
  //   tmp->resize(1,2, false);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(0) == 1, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(1) == 2, true);
  //   for(unsigned int i = 0; i<1; ++i)
  //     for(unsigned int j=0;j<2;++j)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", (*tmp)(i,j) == 0 , true);
  delete tmp;
  cout << "--> resize test ended with success." << endl;
}

void SimpleMatrixTest::testNormInf()
{
  cout << "--> Test: normInf." << endl;
  double n = SicM->normInf();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf: ", n == 7 , true);
  cout << "--> normInf test ended with success." << endl;
}

void SimpleMatrixTest::testSetBlock()
{
  cout << "--> Test: testSetBlock." << endl;

  // Copy of a sub-block of a Simple into a Simple
  SiconosMatrix * MIn = new SimpleMatrix(10, 10);
  for (unsigned int i = 0; i < 10; ++i)
    for (unsigned int j = 0 ; j < 10; ++j)
      (*MIn)(i, j) = i + j;

  SiconosMatrix * MOut = new SimpleMatrix(5, 5);

  Index subDim(2);
  Index subPos(4);
  subDim[0] = 2;
  subDim[1] = 3;
  subPos[0] = 1;
  subPos[1] = 2;
  subPos[2] = 1;
  subPos[3] = 2;

  setBlock(MIn, MOut, subDim, subPos);

  for (unsigned int i = subPos[2]; i < subPos[2] + subDim[0]; ++i)
    for (unsigned int j = subPos[3] ; j < subPos[3] + subDim[1]; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock: ", fabs((*MOut)(i, j) - (*MIn)(i, j)) < tol , true);

  delete MOut;

  // Copy of a sub-block of a Simple into a Block
  Cb->zero();
  setBlock(MIn, Cb, subDim, subPos);

  for (unsigned int i = subPos[2]; i < subPos[2] + subDim[0]; ++i)
    for (unsigned int j = subPos[3] ; j < subPos[3] + subDim[1]; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock: ", fabs((*Cb)(i, j) - (*MIn)(i, j)) < tol , true);
  delete MIn;

  // Copy of a sub-block of a Block into a Simple

  MOut = new SimpleMatrix(5, 5);
  setBlock(Ab, MOut, subDim, subPos);

  for (unsigned int i = subPos[2]; i < subPos[2] + subDim[0]; ++i)
    for (unsigned int j = subPos[3] ; j < subPos[3] + subDim[1]; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock: ", fabs((*MOut)(i, j) - (*Ab)(i, j)) < tol , true);

  delete MOut;
  cout << "-->  setBlock test ended with success." << endl;
}

void SimpleMatrixTest::testSetBlock2()
{
  cout << "--> Test: testSetBlock2." << endl;
  // Copy of a Simple into a sub-block of Simple
  SimpleMatrix * MOut = new SimpleMatrix(10, 10);

  SiconosMatrix * MIn = new SimpleMatrix(5, 5);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0 ; j < 5; ++j)
      (*MIn)(i, j) = i + j;

  MOut->setBlock(2, 3, *MIn);

  for (unsigned int i = 2; i < 7; ++i)
    for (unsigned int j = 3 ; j < 8; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock2: ", fabs((*MOut)(i, j) - (*MIn)(i - 2, j - 3)) < tol , true);

  delete MIn;

  // Copy of a Block into a sub-block of Simple

  MIn = new BlockMatrix(m4, m4, m4, m4);
  MOut->setBlock(2, 3, *MIn);

  for (unsigned int i = 2; i < 6; ++i)
    for (unsigned int j = 3 ; j < 7; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock2: ", fabs((*MOut)(i, j) - (*MIn)(i - 2, j - 3)) < tol , true);

  delete MIn;
  delete MOut;
  cout << "-->  setBlock2 test ended with success." << endl;
}


void SimpleMatrixTest::testGetSetRowCol()
{
  cout << "--> Test: get, set Row and Col." << endl;

  SiconosVector * vIn = new SimpleVector(10, 1.2);
  SiconosVector * vBIn = new BlockVector();
  SiconosVector * v1 = new SimpleVector(3, 2);
  SiconosVector * v2 = new SimpleVector(5, 3);
  SiconosVector * v3 = new SimpleVector(2, 4);
  vBIn->insertPtr(v1);
  vBIn->insertPtr(v2);
  vBIn->insertPtr(v3);

  // Set row with a SimpleVector
  C->setRow(4, *vIn);
  for (unsigned int i = 0; i < C->size(1); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(4, i) - 1.2) < tol, true);

  // Set col with a SimpleVector
  C->setCol(4, *vIn);
  for (unsigned int i = 0; i < C->size(0); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(i, 4) - 1.2) < tol, true);

  // Set row with a BlockVector
  C->setRow(4, *vBIn);

  for (unsigned int i = 0; i < C->size(1); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(4, i) - (*vBIn)(i)) < tol, true);

  // Set col with a BlockVector
  C->setCol(4, *vBIn);
  for (unsigned int i = 0; i < C->size(0); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(i, 4) - (*vBIn)(i)) < tol, true);

  *C = *A; //reset C
  vIn->zero();
  vBIn->zero();
  // get row and copy it into a SimpleVector
  C->getRow(4, *vIn);
  for (unsigned int i = 0; i < C->size(1); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(4, i) - (*vIn)(i)) < tol, true);

  // get col and copy it into a SimpleVector
  C->getCol(4, *vIn);
  for (unsigned int i = 0; i < C->size(0); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(i, 4) - (*vIn)(i)) < tol, true);

  // get row and copy it into a BlockVector
  C->getRow(4, *vBIn);
  for (unsigned int i = 0; i < C->size(1); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(4, i) - (*vBIn)(i)) < tol, true);

  // get col and copy it into a BlockVector
  C->getCol(4, *vBIn);
  for (unsigned int i = 0; i < C->size(0); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(i, 4) - (*vBIn)(i)) < tol, true);

  delete vBIn;
  delete vIn;
  delete v1;
  delete v2;
  delete v3;

  cout << "--> get, set Row and Col tests ended with success." << endl;
}

void SimpleMatrixTest::testTrans()
{
  cout << "--> Test: trans." << endl;

  // Transpose in place ...
  SimpleMatrix * ref = new SimpleMatrix(*D);
  SimpleMatrix * tRef = new SimpleMatrix(*ref);

  tRef->trans();
  for (unsigned int i = 0; i < ref->size(0); ++i)
    for (unsigned int j = 0 ; j < ref->size(1); ++j)
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(j, i) , true);

  // Transpose of another matrix ...
  // Dense
  tRef->zero();

  tRef->trans(*ref);
  for (unsigned int i = 0; i < ref->size(0); ++i)
    for (unsigned int j = 0 ; j < ref->size(1); ++j)
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(j, i) , true);

  delete tRef;
  delete ref;
  // Sym
  ref = new SimpleMatrix(*S);
  tRef = new SimpleMatrix(*ref);
  tRef->trans(*ref);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef) == (*ref) , true);
  delete tRef;
  delete ref;
  // Sparse
  ref = new SimpleMatrix(*SP);
  tRef = new SimpleMatrix(*ref);
  tRef->trans(*ref);
  //   for(unsigned int i = 0; i<ref->size(0); ++i)
  //     {
  //       for(unsigned int j = 0 ; j< ref->size(1); ++j)
  //  if(i==j)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i,j) == (*ref)(i,j) , true);
  //  else
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i,j) == (*ref)(j,i) , true);
  //     }
  delete tRef;
  delete ref;
  // Banded
  //   ref = new SimpleMatrix(*Band);
  //   tRef = new SimpleMatrix(*ref);
  //   *tRef = trans(*ref);
  //   for(unsigned int i = 0; i<ref->size(0); ++i)
  //     for(unsigned int j = 0 ; j< ref->size(1); ++j)
  //       if(i==j)
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i,j) == (*ref)(i,j) , true);
  //       else
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i,j) == (*ref)(j,i) , true);
  //   delete tRef;
  //   delete ref;
  cout << "-->  test trans ended with success." << endl;
}


void SimpleMatrixTest::testAssignment0()
{
  cout << "--> Test: assignment0." << endl;

  // Simple = Simple

  SimpleMatrix * ref = new SimpleMatrix(*D);
  SiconosMatrix * tRef = new SimpleMatrix(*SicM);
  // Dense = any type:
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*T);
  SiconosMatrix * tRef3 = new SimpleMatrix(3, 3);
  *tRef3 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef3) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*S);
  *tRef3 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef3) == (*ref) , true);
  delete ref;

  SiconosMatrix * tRef4 = new SimpleMatrix(4, 4);
  ref = new SimpleMatrix(*SP);
  *tRef4 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef4) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*Band);
  *tRef4 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef4) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*Z);
  *tRef3 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef3) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*I);
  *tRef3 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef3) == (*ref) , true);
  delete ref;

  delete tRef;
  // Triang = Triang, Zero or Identity
  ref = new SimpleMatrix(*T);
  tRef = new SimpleMatrix(*T);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*Z);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*I);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);
  delete ref;

  delete tRef;
  // Sym = Sym, Zero or Id
  ref = new SimpleMatrix(*S);
  tRef = new SimpleMatrix(*S);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);
  delete ref;
  ref = new SimpleMatrix(*Z);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);
  delete ref;
  ref = new SimpleMatrix(*I);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);
  delete ref;
  delete tRef;
  // Sparse = Sparse or Zero
  ref = new SimpleMatrix(*SP);
  tRef = new SimpleMatrix(*SP);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);

  delete ref;
  ref = new SimpleMatrix(*Z2);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);
  delete ref;
  delete tRef;
  // Banded = Banded, Id or Zero
  ref = new SimpleMatrix(*Band);
  tRef = new SimpleMatrix(*Band);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);

  delete ref;
  ref = new SimpleMatrix(*Z2);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);
  delete ref;
  ref = new SimpleMatrix(*I2);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref) , true);
  delete ref;

  delete tRef;
  delete tRef3;
  delete tRef4;

  cout << "-->  test assignment0 ended with success." << endl;
}

void SimpleMatrixTest::testAssignment1()
{
  cout << "--> Test: assignment1." << endl;

  // Simple = Siconos(Block)

  *C = *Ab;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment1: ", (*C) == (*Ab) , true);
  cout << "-->  test assignment1 ended with success." << endl;
}

void SimpleMatrixTest::testAssignment2()
{
  cout << "--> Test: assignment2." << endl;

  // Simple = Siconos(Simple)

  SiconosMatrix * ref = new SimpleMatrix(*D);
  SiconosMatrix * tRef = new SimpleMatrix(*SicM);
  // Dense = any type:
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*T);
  SiconosMatrix * tRef3 = new SimpleMatrix(3, 3);
  *tRef3 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef3) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*S);
  *tRef3 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef3) == (*ref) , true);
  delete ref;

  SiconosMatrix * tRef4 = new SimpleMatrix(4, 4);
  ref = new SimpleMatrix(*SP);
  *tRef4 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef4) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*Band);
  *tRef4 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef4) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*Z);
  *tRef3 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef3) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*I);
  *tRef3 = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef3) == (*ref) , true);
  delete ref;

  delete tRef;
  // Triang = Triang, Zero or Identity
  ref = new SimpleMatrix(*T);
  tRef = new SimpleMatrix(*T);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*Z);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);
  delete ref;

  ref = new SimpleMatrix(*I);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);
  delete ref;

  delete tRef;
  // Sym = Sym, Zero or Id
  ref = new SimpleMatrix(*S);
  tRef = new SimpleMatrix(*S);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);
  delete ref;
  ref = new SimpleMatrix(*Z);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);
  delete ref;
  ref = new SimpleMatrix(*I);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);
  delete ref;
  delete tRef;
  // Sparse = Sparse or Zero
  ref = new SimpleMatrix(*SP);
  tRef = new SimpleMatrix(*SP);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);

  delete ref;
  ref = new SimpleMatrix(*Z2);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);
  delete ref;
  delete tRef;
  // Banded = Banded, Id or Zero
  ref = new SimpleMatrix(*Band);
  tRef = new SimpleMatrix(*Band);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);

  delete ref;
  ref = new SimpleMatrix(*Z2);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);
  delete ref;
  ref = new SimpleMatrix(*I2);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref) , true);
  delete ref;

  delete tRef;
  delete tRef3;
  delete tRef4;
  cout << "-->  test assignment2 ended with success." << endl;
}

void SimpleMatrixTest::testOperators1()
{
  cout << "--> Test: operators1." << endl;
  //+=, -=, *=, /=

  SiconosMatrix * tmp = new SimpleMatrix(*D);
  // Dense *=, /=
  double a = 2.2;
  int a1 = 2;
  *tmp *= a;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*tmp)(i, j) - a * (*D)(i, j)) < tol  , true);

  *tmp *= a1;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*tmp)(i, j) - a * a1 * (*D)(i, j)) < tol  , true);

  *tmp /= a;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*tmp)(i, j) - a1 * (*D)(i, j)) < tol  , true);

  *tmp /= a1;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*tmp)(i, j) - (*D)(i, j)) < tol  , true);

  // Dense +=, -= Dense

  *tmp += *SicM;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*tmp)(i, j) - (*SicM)(i, j) - (*D)(i, j)) < tol , true);

  *tmp -= *SicM;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*tmp)(i, j) - (*D)(i, j)) < tol , true);

  delete tmp;

  // Dense +=, -= Block
  C->zero();
  *C += *Ab;
  *C += *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*C)(i, j) - 2 * (*Ab)(i, j)) < tol , true);
  *C -= *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*C)(i, j) - (*Ab)(i, j)) < tol , true);

  cout << "-->  test operators1 ended with success." << endl;
}

void SimpleMatrixTest::testOperators2()
{
  cout << "--> Test: operators2." << endl;
  // +=, -=, *=, /= triangular
  SiconosMatrix * tmp = new SimpleMatrix(*T);
  SiconosMatrix * tmp2 = new SimpleMatrix(*T);
  *tmp += *tmp2;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.0 * (*T)(i, j) , true);

  int mult = 2;
  double mult0 = 2.2;
  *tmp *= mult0;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*T)(i, j)) < tol , true);

  *tmp *= mult;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*T)(i, j)) < tol , true);

  *tmp /= mult;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*T)(i, j)) < tol , true);

  *tmp /= mult0;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*T)(i, j) , true);

  *tmp -= *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getTriang() - *T) == 0 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", tmp->getNum() == 2 , true);

  delete tmp;
  delete tmp2;
  cout << "-->  test operators2 ended with success." << endl;
}

void SimpleMatrixTest::testOperators3()
{
  cout << "--> Test: operators3." << endl;
  // +=, -=, *=, /= Symmetric
  SiconosMatrix * tmp = new SimpleMatrix(*S);
  SiconosMatrix * tmp2 = new SimpleMatrix(*S);
  *tmp += *tmp2;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.0 * (*S)(i, j) , true);

  int mult = 2;
  double mult0 = 2.2;
  *tmp *= mult0;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*S)(i, j)) < tol , true);

  *tmp *= mult;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*S)(i, j)) < tol , true);

  *tmp /= mult;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*S)(i, j)) < tol , true);

  *tmp /= mult0;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*S)(i, j) , true);

  *tmp -= *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getSym() - *S) == 0 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", tmp->getNum() == 3 , true);

  delete tmp;
  delete tmp2;
  cout << "-->  test operators3 ended with success." << endl;
}

void SimpleMatrixTest::testOperators4()
{
  cout << "--> Test: operators4." << endl;
  // +=, -=, *=, /= sparse
  SiconosMatrix * tmp = new SimpleMatrix(*SP);
  SiconosMatrix * tmp2 = new SimpleMatrix(*SP);
  SiconosMatrix * tmp3 = new SimpleMatrix(*T2);

  SiconosMatrix * tmp4 = new SimpleMatrix(*Band);
  SiconosMatrix * tmp5 = new SimpleMatrix(*S2);

  *tmp += *tmp2;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * (*SP)(i, j)) < tol , true);

  int mult = 2;
  double mult0 = 2.2;
  *tmp *= mult0;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*SP)(i, j)) < tol , true);

  *tmp *= mult;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*SP)(i, j)) < tol , true);

  *tmp /= mult;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*SP)(i, j)) < tol , true);

  *tmp /= mult0;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0 ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2 * (*SP)(i, j)) < tol , true);

  *tmp -= *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getSparse() - *SP) == 0 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", tmp->getNum() == 4 , true);

  // += -= a triangular
  *tmp += *tmp3;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
  {
    for (unsigned int j = 0; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP)(i, j)) < tol , true);
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP)(i, j) - (*tmp3)(i, j)) < tol , true);
  }

  *tmp -= *tmp3;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP)(i, j)) < tol , true);

  // += -= a banded
  *tmp -= *tmp;
  *tmp += *tmp4;
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*Band)(i, j)) < tol , true);

  *tmp -= *tmp4;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) < tol , true);

  // += -= a sym

  *tmp += *tmp5;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
  {
    for (unsigned int j = 0; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP)(i, j) - (*tmp5)(j, i)) < tol , true);
    for (unsigned int j = i ; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP)(i, j) - (*tmp5)(i, j)) < tol , true);
  }

  *tmp -= *tmp5;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP)(i, j)) < tol , true);

  delete tmp;
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  cout << "-->  test operators4 ended with success." << endl;
}

void SimpleMatrixTest::testOperators5()
{
  cout << "--> Test: operators5." << endl;
  // +=, -=, *=, /= banded
  SiconosMatrix * tmp = new SimpleMatrix(*Band);
  SiconosMatrix * tmp2 = new SimpleMatrix(*Band);
  *tmp += *tmp2;
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.0 * (*Band)(i, j) , true);

  int mult = 2;
  double mult0 = 2.2;
  *tmp *= mult0;
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*Band)(i, j)) < tol , true);

  *tmp *= mult;
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*Band)(i, j)) < tol , true);

  *tmp /= mult;
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*Band)(i, j)) < tol , true);

  *tmp /= mult0;
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*Band)(i, j) , true);

  *tmp -= *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getBanded() - *Band) == 0 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", tmp->getNum() == 5 , true);

  delete tmp;
  delete tmp2;
  cout << "-->  test operators5 ended with success." << endl;
}

void SimpleMatrixTest::testOperators6()
{
  cout << "--> Test: operator6." << endl;

  // ============= C = A + B =============

  // Dense = Dense + Dense
  C->zero();
  *C = *A + *B;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*C)(i, j) - (*A)(i, j) - (*B)(i, j)) < tol, true);

  C->zero();
  // Dense = Dense + Block
  *C = *A + *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*C)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();
  // Dense = Block + Dense
  *C = *Ab + *A;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*C)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();

  // Dense = Block + Block
  *C = *Ab + *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*C)(i, j) - (*Ab)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();

  // Block = Dense + Dense
  Cb->zero();
  *Cb = *A + *B;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*Cb)(i, j) - (*A)(i, j) - (*B)(i, j)) < tol, true);

  // Block = Dense + Block

  Cb->zero();
  *Cb = *A + *Ab;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*Cb)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);


  // Block = Block + Dense

  Cb->zero();
  *Cb = *Ab + *A;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*A)(i, j)) < tol, true);


  // Block = Block + Block

  Cb->zero();
  *Cb = *Ab + *Bb;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*Bb)(i, j)) < tol, true);
  Cb->zero();

  // ============= C = A - B =============

  // Dense = Dense - Dense
  C->zero();
  *C = *A - *B;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*C)(i, j) - (*A)(i, j) + (*B)(i, j)) < tol, true);

  C->zero();
  // Dense = Dense - Block
  *C = *A - *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*C)(i, j) - (*A)(i, j) + (*Ab)(i, j)) < tol, true);

  C->zero();
  // Dense = Block - Dense
  *C = *Ab - *A;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*C)(i, j) + (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();

  // Dense = Block - Block
  *C = *Ab - *Bb;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*C)(i, j) - (*Ab)(i, j) + (*Bb)(i, j)) < tol, true);

  C->zero();

  // Block = Dense - Dense
  Cb->zero();
  *Cb = *A - *B;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*Cb)(i, j) - (*A)(i, j) + (*B)(i, j)) < tol, true);

  // Block = Dense - Block

  Cb->zero();
  *Cb = *A - *Ab;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*Cb)(i, j) - (*A)(i, j) + (*Ab)(i, j)) < tol, true);


  // Block = Block - Dense

  Cb->zero();
  *Cb = *Ab - *A;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*Cb)(i, j) - (*Ab)(i, j) + (*A)(i, j)) < tol, true);


  // Block = Block - Block

  Cb->zero();
  *Cb = *Ab - *Bb;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6: ", fabs((*Cb)(i, j) - (*Ab)(i, j) + (*Bb)(i, j)) < tol, true);
  Cb->zero();
  cout << "-->  test operators6 ended with success." << endl;
}

void SimpleMatrixTest::testOperators6Bis()
{
  cout << "--> Test: operator6Bis." << endl;

  // ============= C = A + B =============

  // Dense = Dense + Dense
  C->zero();
  add(*A, *B, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*C)(i, j) - (*A)(i, j) - (*B)(i, j)) < tol, true);

  C->zero();
  // Dense = Dense + Block
  add(*A, *Ab, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*C)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();
  // Dense = Block + Dense
  add(*Ab, *A, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*C)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();

  // Dense = Block + Block
  add(*Ab, *Ab, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*C)(i, j) - (*Ab)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();

  // Block = Dense + Dense
  Cb->zero();
  add(*A, *B, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*Cb)(i, j) - (*A)(i, j) - (*B)(i, j)) < tol, true);

  // Block = Dense + Block

  Cb->zero();
  add(*A, *Ab, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*Cb)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);


  // Block = Block + Dense

  Cb->zero();
  add(*Ab, *A, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*A)(i, j)) < tol, true);


  // Block = Block + Block

  Cb->zero();
  add(*Ab, *Bb, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*Bb)(i, j)) < tol, true);
  Cb->zero();

  // ============= C = A - B =============

  // Dense = Dense - Dense
  C->zero();
  sub(*A, *B, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*C)(i, j) - (*A)(i, j) + (*B)(i, j)) < tol, true);

  C->zero();
  // Dense = Dense - Block
  sub(*A, *Ab, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*C)(i, j) - (*A)(i, j) + (*Ab)(i, j)) < tol, true);

  C->zero();
  // Dense = Block - Dense
  sub(*Ab, *A, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*C)(i, j) + (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();

  // Dense = Block - Block
  sub(*Ab, *Bb, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0 ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*C)(i, j) - (*Ab)(i, j) + (*Bb)(i, j)) < tol, true);

  C->zero();

  // Block = Dense - Dense
  Cb->zero();
  sub(*A, *B, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*Cb)(i, j) - (*A)(i, j) + (*B)(i, j)) < tol, true);

  // Block = Dense - Block

  Cb->zero();
  sub(*A, *Ab, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*Cb)(i, j) - (*A)(i, j) + (*Ab)(i, j)) < tol, true);


  // Block = Block - Dense

  Cb->zero();
  sub(*Ab, *A, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*Cb)(i, j) - (*Ab)(i, j) + (*A)(i, j)) < tol, true);


  // Block = Block - Block

  Cb->zero();
  sub(*Ab, *Bb, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0 ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Bis: ", fabs((*Cb)(i, j) - (*Ab)(i, j) + (*Bb)(i, j)) < tol, true);
  Cb->zero();
  cout << "-->  test operators6Bis ended with success." << endl;
}

void SimpleMatrixTest::testOperators6Ter()
{
  cout << "--> Test: operator6Ter." << endl;

  // +, - , prod for non-dense matrices.

  // Triang +,-,* Triang
  SiconosMatrix * tmp = new SimpleMatrix(*T);
  SiconosMatrix * tmp2 = new SimpleMatrix(*T);
  SiconosMatrix * res = new SimpleMatrix(3, 3, TRIANGULAR);
  *res = *tmp + *tmp2;
  for (unsigned int i = 0; i < res->size(0); ++i)
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res)(i, j) == ((*T)(i, j) + (*T)(i, j)), true);

  *res = *tmp - *tmp2;
  for (unsigned int i = 0; i < res->size(0); ++i)
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res)(i, j) == 0, true);

  *res = prod(*tmp, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", norm_inf(res->getTriang() - prod(*T, *T)) == 0, true);

  delete tmp;
  delete tmp2;
  delete res;

  // Sym +,-,* Sym
  tmp = new SimpleMatrix(*S);
  tmp2 = new SimpleMatrix(*S);
  res = new SimpleMatrix(3, 3, SYMMETRIC);
  *res = *tmp + *tmp2;
  for (unsigned int i = 0; i < res->size(0); ++i)
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res)(i, j) == ((*S)(i, j) + (*S)(i, j)), true);

  *res = *tmp - *tmp2;
  for (unsigned int i = 0; i < res->size(0); ++i)
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res)(i, j) == 0, true);

  *res = prod(*tmp , *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", norm_inf(res->getSym() - prod(*S, *S)) == 0, true);

  delete tmp;
  delete tmp2;
  delete res;

  // Sparse +,-,* Sparse
  tmp = new SimpleMatrix(*SP);
  tmp2 = new SimpleMatrix(*SP);
  res = new SimpleMatrix(4, 4, SPARSE);
  *res = *tmp + *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res) == (2.0 * (*tmp)), true);

  *res = prod(*tmp , *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", norm_inf(*res->getSparsePtr() - prod(*SP, *SP)) < tol, true);

  *res = *tmp - *tmp2;
  tmp->zero();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res) == *tmp , true);

  delete tmp;
  delete tmp2;
  delete res;

  // Banded +,- Banded
  tmp = new SimpleMatrix(*Band);
  tmp2 = new SimpleMatrix(*Band);
  res = new SimpleMatrix(4, 4, BANDED);
  *res = *tmp + *tmp2;
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res)(i, j) == ((*Band)(i, j) + (*Band)(i, j)), true);
  *res = *tmp - *tmp2;
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res)(i, j) == 0, true);

  delete tmp;
  delete tmp2;
  delete res;

  cout << "-->  test operators6Ter6 ended with success." << endl;
}

void SimpleMatrixTest::testOperators7()
{
  cout << "--> Test: operator7." << endl;
  SiconosMatrix * tmp1 = new SimpleMatrix(*D);
  tmp1->resize(4, 4);
  SiconosMatrix * tmp2 = new SimpleMatrix(*T2);
  SiconosMatrix * tmp3 = new SimpleMatrix(*S2);
  SiconosMatrix * tmp4 = new SimpleMatrix(*SP);
  SiconosMatrix * tmp5 = new SimpleMatrix(*Band);
  SiconosMatrix * tmp6 = new SimpleMatrix(*Z2);
  SiconosMatrix * tmp7 = new SimpleMatrix(*I2);

  SiconosMatrix * res = new SimpleMatrix(4, 4);

  // dense + ...
  // ... triang
  add(*tmp1, * tmp2, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp2)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol, true);
  }
  // ... Sym
  add(*tmp1, * tmp3, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp3)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp3)(j, i)) < tol, true);
  }
  // ... Sparse
  add(*tmp1, * tmp4, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
    for (unsigned int j = 0 ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp4)(i, j)) < tol, true);
  // ... Banded
  add(*tmp1, * tmp5, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*Band)(i, j)) < tol, true);
  // Zero
  add(*tmp1, * tmp6, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp1, true);

  // Id
  add(*tmp1, * tmp7, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = 0 ; j < res->size(1); ++j)
    {
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - 1) < tol, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol, true);
    }
  }

  // dense - ...
  // ... triangular
  sub(*tmp1, * tmp2, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + (*tmp2)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol, true);
  }
  // ... Sym
  sub(*tmp1, * tmp3, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + (*tmp3)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + (*tmp3)(j, i)) < tol, true);
  }
  // ... Sparse
  sub(*tmp1, * tmp4, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
    for (unsigned int j = 0 ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + (*tmp4)(i, j)) < tol, true);
  // ... Banded
  sub(*tmp1, * tmp5, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + (*Band)(i, j)) < tol, true);

  // Zero
  sub(*tmp1, * tmp6, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp1, true);

  // Id
  sub(*tmp1, * tmp7, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = 0 ; j < res->size(1); ++j)
    {
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + 1) < tol, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol, true);
    }
  }
  // triang + ...
  // ... dense
  add(*tmp2, * tmp1, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp2)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol, true);
  }
  // ... Sym
  add(*tmp2, * tmp3, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) - (*tmp3)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i)) < tol, true);
  }
  // ... Sparse
  add(*tmp2, * tmp4, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp2)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j)) < tol, true);
  }

  // ... Banded
  add(*tmp2, * tmp5, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) - (*Band)(i, j)) < tol, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*Band)(i, j)) < tol, true);

  // ... Zero
  add(*tmp2, * tmp6, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp2, true);

  // ... Identity
  add(*tmp2, * tmp7, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) - (*tmp7)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j)) < tol, true);
  }

  // triang - ...
  // ... dense
  sub(*tmp2, * tmp1, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) + (*tmp1)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp1)(i, j)) < tol, true);
  }
  // ... Sym
  sub(*tmp2, * tmp3, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) + (*tmp3)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp3)(j, i)) < tol, true);
  }
  // ... Sparse
  sub(*tmp2, * tmp4, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) + (*tmp4)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp4)(i, j)) < tol, true);
  }

  // ... Banded
  sub(*tmp2, * tmp5, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) + (*Band)(i, j)) < tol, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*Band)(i, j)) < tol, true);

  // ... Zero
  sub(*tmp2, * tmp6, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp2, true);

  // Identity
  sub(*tmp2, * tmp7, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) + (*tmp7)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j)) < tol, true);
  }

  // sym + ...
  // ... dense
  add(*tmp3, * tmp1, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp3)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp3)(j, i)) < tol, true);
  }
  // ... triang
  add(*tmp3, * tmp2, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) - (*tmp2)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i)) < tol, true);
  }
  // ... Sparse
  add(*tmp3, * tmp4, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp3)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp3)(j, i)) < tol, true);
  }

  // ... Banded
  add(*tmp3, * tmp5, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) - (*Band)(i, j)) < tol, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) - (*Band)(i, j)) < tol, true);


  // ... Zero
  add(*tmp3, * tmp6, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp3, true);

  // ... identity
  add(*tmp3, * tmp7, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp7)(i, j) - (*tmp3)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp7)(i, j) - (*tmp3)(j, i)) < tol, true);
  }

  // sym - ...
  // ... dense
  sub(*tmp3, * tmp1, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) + (*tmp1)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) + (*tmp1)(i, j)) < tol, true);
  }
  // ... triang
  sub(*tmp3, * tmp2, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) + (*tmp2)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i)) < tol, true);
  }
  // ... Sparse
  sub(*tmp3, * tmp4, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) + (*tmp4)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) + (*tmp4)(i, j)) < tol, true);
  }

  // ... Banded
  sub(*tmp3, * tmp5, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) + (*Band)(i, j)) < tol, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) + (*Band)(i, j)) < tol, true);

  // ... Zero
  sub(*tmp3, * tmp6, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: "  , *res == *tmp3, true);
  // Identity
  sub(*tmp3, * tmp7, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) + (*tmp7)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) + (*tmp7)(i, j)) < tol, true);
  }

  // sparse + ...
  // ... dense
  add(*tmp4, * tmp1, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
    for (unsigned int j = 0 ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp4)(i, j)) < tol, true);
  // ... triang
  add(*tmp4, * tmp2, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp2)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j)) < tol, true);
  }
  // ... Sym
  add(*tmp4, * tmp3, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp3)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp3)(j, i)) < tol, true);
  }
  // ... Banded
  add(*tmp4, * tmp5, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*Band)(i, j)) < tol, true);

  // ... zero
  add(*tmp4, * tmp6, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: "  , *res == *tmp4, true);

  // sparse - ...
  // ... dense
  sub(*tmp4, * tmp1, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
    for (unsigned int j = 0 ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) + (*tmp1)(i, j)) < tol, true);
  // ... triangular
  sub(*tmp4, * tmp2, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) + (*tmp2)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j)) < tol, true);
  }
  // ... Sym
  sub(*tmp4, * tmp3, *res);
  for (unsigned int i = 0; i < res->size(0); ++i)
  {
    for (unsigned int j = i ; j < res->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) + (*tmp3)(i, j)) < tol, true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) + (*tmp3)(j, i)) < tol, true);
  }

  // ... Banded
  sub(*tmp4, * tmp5, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) + (*Band)(i, j)) < tol, true);

  // ... zero
  sub(*tmp4, * tmp6, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: "  , *res == *tmp4, true);

  // Banded + ...
  // ... dense
  add(*tmp5, * tmp1, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol, true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp5)(i, j)) < tol, true);
    for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol, true);
  }
  // ... triang
  add(*tmp5, * tmp2, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j)) < tol, true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) - (*tmp5)(i, j)) < tol, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j)) < tol, true);
    for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j)) < tol, true);
  }

  // ...sym
  add(*tmp5, * tmp3, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j)) < tol, true);
      else  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i)) < tol, true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) - (*tmp5)(i, j)) < tol, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) - (*tmp5)(i, j)) < tol, true);
    for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j)) < tol, true);
      else  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i)) < tol, true);
  }

  //... sparse
  add(*tmp5, * tmp4, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j)) < tol, true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp5)(i, j)) < tol, true);
    for (signed j = std::min(i + 2, signed(Band->size1())); j < signed(Band->size1()); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j)) < tol, true);
  }

  // ... zero
  add(*tmp5, * tmp6, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: "  , *res == *tmp5, true);
  // ... identity
  add(*tmp5, * tmp7, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp7)(i, j)) < tol, true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp7)(i, j) - (*tmp5)(i, j)) < tol, true);
    for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp7)(i, j)) < tol, true);
  }

  // Banded - ...
  // ... dense

  sub(*tmp5, * tmp1, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp1)(i, j)) < tol, true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j) + (*tmp1)(i, j)) < tol, true);
    for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp1)(i, j)) < tol, true);
  }

  // ... triang
  sub(*tmp5, * tmp2, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ",  fabs((*res)(i, j) + (*tmp2)(i, j)) < tol, true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j) + (*tmp2)(i, j)) < tol, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j)) < tol, true);
    for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp2)(i, j)) < tol, true);
  }

  // ...sym
  sub(*tmp5, * tmp3, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp3)(i, j)) < tol, true);
      else  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp3)(j, i)) < tol, true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j) + (*tmp3)(i, j)) < tol, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j) + (*tmp3)(j, i)) < tol, true);
    for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp3)(i, j)) < tol, true);
      else  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp3)(j, i)) < tol, true);
  }

  //... sparse
  sub(*tmp5, * tmp4, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp4)(i, j)) < tol, true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp4)(i, j) - (*tmp5)(i, j)) < tol, true);
    for (signed j = std::min(i + 2, signed(Band->size1())); j < signed(Band->size1()); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp4)(i, j)) < tol, true);
  }

  // ... zero
  sub(*tmp5, * tmp6, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: "  , *res == *tmp5, true);
  // ... identity
  sub(*tmp5, * tmp7, *res);
  for (signed i = 0; i < signed(Band->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp7)(i, j)) < tol, true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j) + (*tmp7)(i, j)) < tol, true);
    for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp7)(i, j)) < tol, true);
  }

  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  delete tmp6;
  delete tmp7;
  delete res;

  cout << "-->  test operators7 ended with success." << endl;
}



void SimpleMatrixTest::testOperators8()
{
  cout << "--> Test: operator8." << endl;

  // Simple = Simple * Simple
  *C = prod(*A, *B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(*C->getDensePtr() - prod(*A->getDensePtr(), *B->getDensePtr())) < tol, true);

  // Block = Simple * Simple
  *Cb = prod(*A, *B);
  DenseMat Dtmp = prod(*A->getDensePtr(), *B->getDensePtr());
  SimpleMatrix * tmp = new SimpleMatrix(Dtmp);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*Cb)(i, j) - (*tmp)(i, j)) < tol, true);

  // Others ...

  SiconosMatrix * tmp1 = new SimpleMatrix(4, 4, 2.3);
  SiconosMatrix * tmp2 = new SimpleMatrix(*T2);
  SiconosMatrix * tmp3 = new SimpleMatrix(*S2);
  SiconosMatrix * tmp4 = new SimpleMatrix(*SP);
  SiconosMatrix * tmp5 = new SimpleMatrix(*Band);

  SiconosMatrix * res = new SimpleMatrix(4, 4);

  // Dense * ...
  // triang
  *res = prod(*tmp1, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp2->getTriang())) < tol, true);
  // Sym
  *res = prod(*tmp1, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp3->getSym())) < tol, true);
  // Sparse
  *res = prod(*tmp1, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp4->getSparse())) < tol, true);
  // Banded
  *res = prod(*tmp1, *tmp5);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp5->getBanded())) < tol, true);
  // triang * ...
  // dense
  *res = prod(*tmp2, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp1->getDense())) < tol, true);
  // Sym
  *res = prod(*tmp2, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp3->getSym())) < tol, true);
  // Sparse
  *res = prod(*tmp2, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp4->getSparse())) < tol, true);
  // Banded
  *res = prod(*tmp2, *tmp5);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp5->getBanded())) < tol, true);
  // sym * ...
  // dense
  *res = prod(*tmp3, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp1->getDense())) < tol, true);
  // triang
  *res = prod(*tmp3, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp2->getTriang())) < tol, true);
  // Sparse
  *res = prod(*tmp3, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp4->getSparse())) < tol, true);
  // Banded
  *res = prod(*tmp3, *tmp5);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp5->getBanded())) < tol, true);
  // Sparse * ...
  // dense
  *res = prod(*tmp4, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp1->getDense())) < tol, true);
  // triang
  *res = prod(*tmp4, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp2->getTriang())) < tol, true);
  // Sym
  *res = prod(*tmp4, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp3->getSym())) < tol, true);
  // Banded
  *res = prod(*tmp4, *tmp5);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp5->getBanded())) < tol, true);
  // Banded * ...
  // dense
  *res = prod(*tmp5, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp1->getDense())) < tol, true);
  // triang
  *res = prod(*tmp5, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp2->getTriang())) < tol, true);
  // Sparse
  *res = prod(*tmp5, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp4->getSparse())) < tol, true);
  // Sym
  *res = prod(*tmp5, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp3->getSym())) < tol, true);

  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  delete res;
  cout << "-->  test operators8 ended with success." << endl;
}

void SimpleMatrixTest::testOperators8Bis()
{
  cout << "--> Test: operator8Bis." << endl;
  // Simple = Simple * Simple
  prod(*A, *B, *C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(*C->getDensePtr() - prod(*A->getDensePtr(), *B->getDensePtr())) < tol, true);

  // Block = Simple * Simple
  prod(*A, *B, *Cb);
  DenseMat Dtmp = prod(*A->getDensePtr(), *B->getDensePtr());
  SimpleMatrix * tmp = new SimpleMatrix(Dtmp);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*Cb)(i, j) - (*tmp)(i, j)) < tol, true);

  // Others ...

  // Others ...
  SiconosMatrix * tmp1 = new SimpleMatrix(4, 4, 2.4);
  SiconosMatrix * tmp2 = new SimpleMatrix(*T2);
  SiconosMatrix * tmp3 = new SimpleMatrix(*S2);
  SiconosMatrix * tmp4 = new SimpleMatrix(*SP);
  SiconosMatrix * tmp5 = new SimpleMatrix(*Band);

  SiconosMatrix * res = new SimpleMatrix(4, 4);

  // Dense * ...
  // triang
  prod(*tmp1, *tmp2, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp2->getTriang())) < tol, true);
  // Sym
  prod(*tmp1, *tmp3, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp3->getSym())) < tol, true);
  // Sparse
  prod(*tmp1, *tmp4, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp4->getSparse())) < tol, true);
  // Banded
  prod(*tmp1, *tmp5, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp5->getBanded())) < tol, true);
  // triang * ...
  // dense
  prod(*tmp2, *tmp1, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp1->getDense())) < tol, true);
  // Sym
  prod(*tmp2, *tmp3, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp3->getSym())) < tol, true);
  // Sparse
  prod(*tmp2, *tmp4, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp4->getSparse())) < tol, true);
  // Banded
  prod(*tmp2, *tmp5, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp5->getBanded())) < tol, true);
  // sym * ...
  // dense
  prod(*tmp3, *tmp1, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp1->getDense())) < tol, true);
  // triang
  prod(*tmp3, *tmp2, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp2->getTriang())) < tol, true);
  // Sparse
  prod(*tmp3, *tmp4, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp4->getSparse())) < tol, true);
  // Banded
  prod(*tmp3, *tmp5, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp5->getBanded())) < tol, true);
  // Sparse * ...
  // dense
  prod(*tmp4, *tmp1, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp1->getDense())) < tol, true);
  // triang
  prod(*tmp4, *tmp2, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp2->getTriang())) < tol, true);
  // Sym
  prod(*tmp4, *tmp3, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp3->getSym())) < tol, true);
  // Banded
  prod(*tmp4, *tmp5, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp5->getBanded())) < tol, true);
  // Banded * ...
  // dense
  prod(*tmp5, *tmp1, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp1->getDense())) < tol, true);
  // triang
  prod(*tmp5, *tmp2, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp2->getTriang())) < tol, true);
  // Sparse
  prod(*tmp5, *tmp4, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp4->getSparse())) < tol, true);
  // Sym
  prod(*tmp5, *tmp3, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Bis: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp3->getSym())) < tol, true);

  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  delete res;
  cout << "-->  test operators8Bis ended with success." << endl;
}

void SimpleMatrixTest::testOperators8Ter()
{
  cout << "--> Test: operator8Ter." << endl;
  // Simple = Simple * Simple
  axpy_prod(*A, *B, *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Ter: ", norm_inf(*C->getDensePtr() - prod(*A->getDensePtr(), *B->getDensePtr())) < tol, true);

  // Simple += Simple * Simple
  SiconosMatrix * backUp = new SimpleMatrix(*C);

  axpy_prod(*A, *B, *C, false);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8Ter: ", norm_inf(*C->getDensePtr() - prod(*A->getDensePtr(), *B->getDensePtr()) - *backUp->getDensePtr()) < tol, true);
  // Block = Simple * Simple
  axpy_prod(*A, *B, *Cb, true);
  DenseMat Dtmp = prod(*A->getDensePtr(), *B->getDensePtr());
  SimpleMatrix * tmp = new SimpleMatrix(Dtmp);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*Cb)(i, j) - (*tmp)(i, j)) < tol, true);

  *backUp = *Cb;
  // Block += Simple * Simple
  axpy_prod(*A, *B, *Cb, false);
  Dtmp = prod(*A->getDensePtr(), *B->getDensePtr());
  *tmp = Dtmp;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*Cb)(i, j) - (*tmp)(i, j) - (*backUp)(i, j)) < tol, true);

  delete backUp;
  cout << "-->  test operators8Ter ended with success." << endl;
}

void SimpleMatrixTest::testOperators8_4() // C += A*B
{
  cout << "--> Test: operator8_4." << endl;
  // Simple = Simple * Simple
  C->zero();
  prod(*A, *B, *C, false);
  prod(*A, *B, *C, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_4: ", norm_inf(*C->getDensePtr() - 2 * prod(*A->getDensePtr(), *B->getDensePtr())) < tol, true);

  // Block = Simple * Simple
  Cb->zero();
  prod(*A, *B, *Cb, false);
  prod(*A, *B, *Cb, false);
  DenseMat Dtmp = prod(*A->getDensePtr(), *B->getDensePtr());
  SimpleMatrix * tmp = new SimpleMatrix(Dtmp);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*Cb)(i, j) - 2 * (*tmp)(i, j)) < tol, true);
  cout << "-->  test operators8_4 ended with success." << endl;
}

void SimpleMatrixTest::testOperators8_5()
{
  // == Test subprod ==

  cout << "--> Test: operator8_5." << endl;
  std::vector<unsigned int> coord(8);
  SiconosVector * x1 = new SimpleVector(2);
  SiconosVector * x2 = new SimpleVector(3);
  SiconosVector * x3 = new SimpleVector(5);
  SiconosVector * y = new SimpleVector(size);
  SiconosVector * x = new BlockVector();
  SiconosVector * v = new SimpleVector(size);
  x->insertPtr(x1);
  x->insertPtr(x2);
  x->insertPtr(x3);
  for (unsigned int i = 0 ; i < size; ++i)
  {
    (*x)(i) = (double)i + 3;
    (*v)(i) = (double)i + 3;
  }

  // v == x but x is a 3-blocks vector.

  // Simple = Simple * Simple, all dense
  // subprod but with full matrix/vectors
  coord[0] = 0;
  coord[1] = size;
  coord[2] = 0;
  coord[3] = size;
  coord[4] = 0;
  coord[5] = size;
  coord[6] = 0;
  coord[7] = size;
  subprod(*A, *v, *y, coord, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", norm_inf(*y->getDensePtr() - prod(*A->getDensePtr(), *v->getDensePtr())) < tol, true);

  // Simple = Simple * Block, all dense
  // subprod but with full matrix/vectors
  subprod(*A, *x, *y, coord, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", norm_inf(*y->getDensePtr() - prod(*A->getDensePtr(), *v->getDensePtr())) < tol, true);

  coord[0] = 0;
  coord[1] = 2;
  coord[2] = 1;
  coord[3] = 3;
  coord[4] = 3;
  coord[5] = 5;
  coord[6] = 2;
  coord[7] = 4;
  y->zero();
  // Simple = Simple * Simple, all dense
  subprod(*A, *v, *y, coord, true);
  double res = (*A)(0, 1) * (*v)(3) + (*A)(0, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A)(1, 1) * (*v)(3) + (*A)(1, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  }
  y->zero();
  // Simple = Simple * Block, all dense
  subprod(*A, *x, *y, coord, true);
  res = (*A)(0, 1) * (*x)(3) + (*A)(0, 2) * (*x)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A)(1, 1) * (*x)(3) + (*A)(1, 2) * (*x)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  }
  //   // Others ...
  // Triang

  SiconosMatrix * A2 = new SimpleMatrix(10, 10, TRIANGULAR);
  for (unsigned i = 0; i < A2->size(0); ++ i)
    for (unsigned j = i; j < A2->size(1); ++ j)
      (*A2)(i, j) = 3 * i + j;

  subprod(*A2, *v, *y, coord, true);
  res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  }
  // Sym
  delete A2;
  A2 = new SimpleMatrix(10, 10, SYMMETRIC);
  for (unsigned i = 0; i < A2->size(0); ++ i)
    for (unsigned j = i; j < A2->size(1); ++ j)
      (*A2)(i, j) = 3 * i + j;

  subprod(*A2, *v, *y, coord, true);
  res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  }

  // Sparse
  delete A2;
  A2 = new SimpleMatrix(10, 10, SPARSE);
  for (unsigned i = 0; i < A2->size(0); ++ i)
    for (unsigned j = i; j < A2->size(1); ++ j)
      A2->setValue(i, j, 3 * i + j);

  subprod(*A2, *v, *y, coord, true);
  res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  }

  // Banded
  A2 = new SimpleMatrix(10, 10, BANDED);
  for (signed i = 0; i < signed(A2->size(0)); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(A2->size(1))); ++ j)
      (*A2)(i, j) = 3 * i + j;
  subprod(*A2, *v, *y, coord, true);
  res = (*A2)(0, 1) * (*v)(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  }

  delete A2;
  delete x1;
  delete x2;
  delete x3;
  delete y;
  delete v;
  delete x;
  cout << "-->  test operators8_5 ended with success." << endl;
}

void SimpleMatrixTest::testOperators8_6()
{
  // == Test subprod, with += ==

  cout << "--> Test: operator8_6." << endl;
  std::vector<unsigned int> coord(8);
  SiconosVector * x1 = new SimpleVector(2);
  SiconosVector * x2 = new SimpleVector(3);
  SiconosVector * x3 = new SimpleVector(5);
  SiconosVector * y = new SimpleVector(size);
  SiconosVector * x = new BlockVector();
  SiconosVector * v = new SimpleVector(size);
  x->insertPtr(x1);
  x->insertPtr(x2);
  x->insertPtr(x3);
  for (unsigned int i = 0 ; i < size; ++i)
  {
    (*x)(i) = (double)i + 3;
    (*v)(i) = (double)i + 3;
  }

  // v == x but x is a 3-blocks vector.

  *y = *v;

  // Simple = Simple * Simple, all dense
  // subprod but with full matrix/vectors
  coord[0] = 0;
  coord[1] = size;
  coord[2] = 0;
  coord[3] = size;
  coord[4] = 0;
  coord[5] = size;
  coord[6] = 0;
  coord[7] = size;
  subprod(*A, *v, *y, coord, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", norm_inf(*y->getDensePtr() - prod(*A->getDensePtr(), *v->getDensePtr()) - *v->getDensePtr()) < tol, true);

  // Simple = Simple * Block, all dense
  // subprod but with full matrix/vectors
  *y = *v;
  subprod(*A, *x, *y, coord, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", norm_inf(*y->getDensePtr() - prod(*A->getDensePtr(), *v->getDensePtr()) - *v->getDensePtr()) < tol, true);

  coord[0] = 0;
  coord[1] = 2;
  coord[2] = 1;
  coord[3] = 3;
  coord[4] = 3;
  coord[5] = 5;
  coord[6] = 2;
  coord[7] = 4;

  // Simple = Simple * Simple, all dense
  *y = *v;
  subprod(*A, *v, *y, coord, false);
  double res = (*A)(0, 1) * (*v)(3) + (*A)(0, 2) * (*v)(4) + (*v)(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A)(1, 1) * (*v)(3) + (*A)(1, 2) * (*v)(4) + (*v)(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  }
  *y = *v;
  // Simple = Simple * Block, all dense
  subprod(*A, *x, *y, coord, false);
  res = (*A)(0, 1) * (*x)(3) + (*A)(0, 2) * (*x)(4) + (*v)(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A)(1, 1) * (*x)(3) + (*A)(1, 2) * (*x)(4) + (*v)(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  }

  //   // Others ...
  // Triang

  SiconosMatrix * A2 = new SimpleMatrix(10, 10, TRIANGULAR);
  for (unsigned i = 0; i < A2->size(0); ++ i)
    for (unsigned j = i; j < A2->size(1); ++ j)
      (*A2)(i, j) = 3 * i + j;

  *y = *v;
  subprod(*A2, *v, *y, coord, false);
  res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4) + (*v)(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4) + (*v)(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  }

  // Sym
  delete A2;
  A2 = new SimpleMatrix(10, 10, SYMMETRIC);
  for (unsigned i = 0; i < A2->size(0); ++ i)
    for (unsigned j = i; j < A2->size(1); ++ j)
      (*A2)(i, j) = 3 * i + j;

  *y = *v;
  subprod(*A2, *v, *y, coord, false);
  res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4) + (*v)(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4) + (*v)(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  }

  // Sparse
  delete A2;
  A2 = new SimpleMatrix(10, 10, SPARSE);
  for (unsigned i = 0; i < A2->size(0); ++ i)
    for (unsigned j = i; j < A2->size(1); ++ j)
      A2->setValue(i, j, 3 * i + j);

  *y = *v;
  subprod(*A2, *v, *y, coord, false);
  res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4) + (*v)(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4) + (*v)(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  }

  // Banded
  A2 = new SimpleMatrix(10, 10, BANDED);
  for (signed i = 0; i < signed(A2->size(0)); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(A2->size(1))); ++ j)
      (*A2)(i, j) = 3 * i + j;
  *y = *v;
  subprod(*A2, *v, *y, coord, false);
  res = (*A2)(0, 1) * (*v)(3) + (*v)(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4) + (*v)(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  }

  delete A2;
  delete x1;
  delete x2;
  delete x3;
  delete y;
  delete v;
  delete x;
  cout << "-->  test operators8_6 ended with success." << endl;
}

void SimpleMatrixTest::testOperators9()
{
  cout << "--> Test: operator9." << endl;

  // C = a*A or A/a

  double a = 2.2;
  int a1 = 3;

  // Simple = a * Simple or Simple/a
  *C = a * *A;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - a * (*A)(i, j)) < tol, true);
  *C = *A * a;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - a * (*A)(i, j)) < tol, true);
  *C = *A * a1;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - a1 * (*A)(i, j)) < tol, true);
  *C = a1 * *A;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - a1 * (*A)(i, j)) < tol, true);

  *C = *A / a;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*A)(i, j) / a) < tol, true);
  *C = *A / a1;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*A)(i, j) / a1) < tol, true);

  // Simple = a * Block

  *C = a * *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - a * (*Ab)(i, j)) < tol, true);
  *C = *Ab * a;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - a * (*Ab)(i, j)) < tol, true);
  *C = *Ab * a1;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - a1 * (*Ab)(i, j)) < tol, true);
  *C = a1 * *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - a1 * (*Ab)(i, j)) < tol, true);

  *C = *Ab / a;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*Ab)(i, j) / a) < tol, true);
  *C = *Ab / a1;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*Ab)(i, j) / a1) < tol, true);

  // Block = a * Block
  *Cb = a * *Ab;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - a * (*Ab)(i, j)) < tol, true);
  *Cb = *Ab * a;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - a * (*Ab)(i, j)) < tol, true);
  *Cb = *Ab * a1;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - a1 * (*Ab)(i, j)) < tol, true);
  *Cb = a1 * *Ab;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - a1 * (*Ab)(i, j)) < tol, true);

  *Cb = *Ab / a;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*Ab)(i, j) / a) < tol, true);
  *Cb = *Ab / a1;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*Ab)(i, j) / a1) < tol, true);

  // Block = a * Simple
  *Cb = a * *A;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - a * (*A)(i, j)) < tol, true);
  *Cb = *A * a;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - a * (*A)(i, j)) < tol, true);
  *Cb = *A * a1;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - a1 * (*A)(i, j)) < tol, true);
  *Cb = a1 * *A;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - a1 * (*A)(i, j)) < tol, true);

  *Cb = *A / a;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*A)(i, j) / a) < tol, true);
  *Cb = *A / a1;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*A)(i, j) / a1) < tol, true);
  cout << "-->  test operators9 ended with success." << endl;
}

void SimpleMatrixTest::testOperators9Bis()
{
  cout << "--> Test: operator9Bis." << endl;

  // C = a*A or A/a

  double a = 2.2;

  // Simple = a * Simple or Simple/a
  scal(a, *A, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Bis: ", fabs((*C)(i, j) - a * (*A)(i, j)) < tol, true);

  scal(1.0 / a, *A, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Bis: ", fabs((*C)(i, j) - (*A)(i, j) / a) < tol, true);
  // Simple = a * Block

  scal(a, *Ab, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Bis: ", fabs((*C)(i, j) - a * (*Ab)(i, j)) < tol, true);

  scal(1.0 / a, *Ab, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Bis: ", fabs((*C)(i, j) - (*Ab)(i, j) / a) < tol, true);

  // Block = a * Block
  scal(a, *Ab, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Bis: ", fabs((*Cb)(i, j) - a * (*Ab)(i, j)) < tol, true);

  scal(1.0 / a, *Ab, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Bis: ", fabs((*Cb)(i, j) - (*Ab)(i, j) / a) < tol, true);

  // Block = a * Simple
  scal(a, *A, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Bis: ", fabs((*Cb)(i, j) - a * (*A)(i, j)) < tol, true);

  scal(1.0 / a, *A, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Bis: ", fabs((*Cb)(i, j) - (*A)(i, j) / a) < tol, true);
  cout << "-->  test operators9Bis ended with success." << endl;
}

void SimpleMatrixTest::testOperators9Ter()
{
  cout << "--> Test: operator9Ter." << endl;

  // C += a*A or A/a

  double a = 2.2;
  C->zero();
  // Simple = a * Simple or Simple/a
  scal(a, *A, *C, false);
  scal(a, *A, *C, false);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Ter: ", fabs((*C)(i, j) - 2 * a * (*A)(i, j)) < tol, true);

  // Simple = a * Block
  C->zero();
  scal(a, *Ab, *C, false);
  scal(a, *Ab, *C, false);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i ; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Ter: ", fabs((*C)(i, j) - 2 * a * (*Ab)(i, j)) < tol, true);

  // Block = a * Block
  Cb->zero();
  scal(a, *Ab, *Cb, false);
  scal(a, *Ab, *Cb, false);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Ter: ", fabs((*Cb)(i, j) - 2 * a * (*Ab)(i, j)) < tol, true);

  // Block = a * Simple
  Cb->zero();
  scal(a, *A, *Cb, false);
  scal(a, *A, *Cb, false);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i ; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9Ter: ", fabs((*Cb)(i, j) - 2 * a * (*A)(i, j)) < tol, true);

  cout << "-->  test operators9Ter ended with success." << endl;
}

void SimpleMatrixTest::testOperators10()
{
  cout << "--> Test: operator10." << endl;
  double m = 2.2;
  int i = 3;
  SiconosMatrix * tmp1 = new SimpleMatrix(*T);
  SiconosMatrix * res = new SimpleMatrix(3, 3, TRIANGULAR);
  *res = m ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang()*m) < tol, true);
  *res = i ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang()*i) < tol, true);
  *res = *tmp1 * m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang()*m) < tol, true);
  *res = *tmp1 * i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang()*i) < tol, true);
  *res = *tmp1 / m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang() / m) < tol, true);
  *res = *tmp1 / i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang() / i) < tol, true);
  delete tmp1;
  delete res;
  cout << "-->  test operators10 ended with success." << endl;
}

void SimpleMatrixTest::testOperators11()
{
  cout << "--> Test: operator11." << endl;
  double m = 2.2;
  int i = 3;
  SiconosMatrix * tmp1 = new SimpleMatrix(*S);
  SiconosMatrix * res = new SimpleMatrix(3, 3, SYMMETRIC);
  *res = m ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym()*m) < tol, true);
  *res = i ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym()*i) < tol, true);
  *res = *tmp1 * m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym()*m) < tol, true);
  *res = *tmp1 * i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym()*i) < tol, true);
  *res = *tmp1 / m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym() / m) < tol, true);
  *res = *tmp1 / i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym() / i) < tol, true);
  delete tmp1;
  delete res;
  cout << "-->  test operator11 ended with success." << endl;
}

void SimpleMatrixTest::testOperators12()
{
  cout << "--> Test: operator12." << endl;
  double m = 2.2;
  int i = 3;
  SiconosMatrix * tmp1 = new SimpleMatrix(*SP);
  SiconosMatrix * res = new SimpleMatrix(4, 4, SPARSE);
  *res = m ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse()*m) < tol, true);
  *res = i ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse()*i) < tol, true);
  *res = *tmp1 * m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse()*m) < tol, true);
  *res = *tmp1 * i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse()*i) < tol, true);
  *res = *tmp1 / m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse() / m) < tol, true);
  *res = *tmp1 / i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse() / i) < tol, true);
  delete tmp1;
  delete res;
  cout << "-->  test operators12 ended with success." << endl;
}

void SimpleMatrixTest::testOperators13()
{
  cout << "--> Test: operator13." << endl;
  //   double m = 2.2;
  //   int i = 3;
  //   SiconosMatrix * tmp1 = new SimpleMatrix(*Band);
  //   SiconosMatrix * res = new SimpleMatrix(*Band);//4,4,BANDED,1,1);
  //   *res = m * *tmp1;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()- tmp1->getBanded()*m)<tol, true);
  //   *res = i ** tmp1;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()- tmp1->getBanded()*i)<tol, true);
  //   *res = *tmp1 * m;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()- tmp1->getBanded()*m)<tol, true);
  //   *res = *tmp1 * i;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()- tmp1->getBanded()*i)<tol, true);
  //   *res = *tmp1 / m;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()- tmp1->getBanded()/m)<tol, true);
  //   *res = *tmp1 / i;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()- tmp1->getBanded()/i)<tol, true);
  //   delete tmp1;
  //   delete res;
  cout << "-->  test operators13 ended with success." << endl;
}

void SimpleMatrixTest::testPow()
{
  cout << "--> Test: pow." << endl;
  // Dense
  SiconosMatrix * tmp1 = new SimpleMatrix(*D);
  SiconosMatrix * res = new SimpleMatrix(2, 2);
  *res = pow(*tmp1, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == prod(*tmp1, prod(*tmp1, *tmp1)), true);
  delete res;
  // Triang
  SiconosMatrix * tmp2 = new SimpleMatrix(*T);
  res = new SimpleMatrix(3, 3, TRIANGULAR);
  *res = pow(*tmp2, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == prod(*tmp2, prod(*tmp2, *tmp2)), true);
  delete res;
  // Sym
  SiconosMatrix * tmp3 = new SimpleMatrix(*S);
  res = new SimpleMatrix(3, 3, SYMMETRIC);
  *res = pow(*tmp3, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == prod(*tmp3, prod(*tmp3, *tmp3)), true);
  delete res;
  // Sparse
  SiconosMatrix * tmp4 = new SimpleMatrix(*SP);
  res = new SimpleMatrix(4, 4, SPARSE);
  *res = pow(*tmp4, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == prod(*tmp4, prod(*tmp4, *tmp4)), true);
  delete res;
  // Banded
  SiconosMatrix * tmp5 = new SimpleMatrix(*Band);
  res = new SimpleMatrix(4, 4);
  *res = pow(*tmp5, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == prod(*tmp5, prod(*tmp5, *tmp5)), true);
  delete res;
  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  cout << "-->  test pow ended with success." << endl;
}

void SimpleMatrixTest::testProd() // y = A*x
{
  cout << "--> Test: prod. mat-vect" << endl;

  SiconosVector * y = new SimpleVector(size);
  SiconosVector * x = new SimpleVector(size, 4.3);
  SiconosVector * x1 = new SimpleVector(size - 2, 2.3);
  SiconosVector * x2 = new SimpleVector(2, 3.1);

  SiconosVector * xB = new BlockVector(x1, x2);
  SiconosVector * yB = new BlockVector(*xB);
  yB->zero();

  // Matrix - vector product

  // Simple = Simple * Simple
  *y = prod(*A, *x);
  double sum;
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", fabs((*y)(i) - sum) < tol, true);
  }
  // Simple = Simple * Block
  *y = prod(*A , *xB);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", fabs((*y)(i) - sum) < tol, true);
  }


  // Block = Simple * Simple
  *yB = prod(*A , *x);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", fabs((*yB)(i) - sum) < tol, true);
  }

  // Block = Simple * Block
  *yB = prod(*A , *xB);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", fabs((*yB)(i) - sum) < tol, true);
  }

  delete x;
  delete y;
  delete xB;
  delete yB;
  delete x1;
  delete x2;

  // Others or old stuff ...

  SiconosMatrix * tmp2 = new SimpleMatrix(*T);
  SiconosMatrix * tmp3 = new SimpleMatrix(*S);
  SiconosMatrix * tmp4 = new SimpleMatrix(*SP);
  SiconosMatrix * tmp5 = new SimpleMatrix(*Band2);
  SiconosVector * v = new SimpleVector(3);
  (*v)(0) = 1;
  (*v)(1) = 2;
  (*v)(2) = 3;
  SiconosVector * vv = new SimpleVector(4);
  (*vv)(0) = 1;
  (*vv)(1) = 2;
  (*vv)(2) = 3;
  SparseVect * sv = new SparseVect(3);
  (*sv)(0) = 4;
  (*sv)(1) = 5;
  (*sv)(2) = 6;
  SparseVect * sv2 = new SparseVect(4);
  (*sv2)(0) = 4;
  (*sv2)(1) = 5;
  (*sv2)(2) = 6;
  SiconosVector * w = new SimpleVector(*sv);
  SiconosVector * ww = new SimpleVector(*sv2);
  SiconosVector * res = new SimpleVector(4);
  SiconosVector * res2 = new SimpleVector(3);

  // Triang * ...
  *res2 = prod(*tmp2, *v);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", norm_2(res2->getDense() - prod(tmp2->getTriang(), v->getDense())) < tol, true);
  *res2 = prod(*tmp2, *w);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", norm_2(res2->getDense() - prod(tmp2->getTriang(), w->getSparse())) < tol, true);
  //   Sym * ...
  *res2 = prod(*tmp3, *v);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", norm_2(res2->getDense() - prod(tmp3->getSym(), v->getDense())) < tol, true);
  *res2 = prod(*tmp3, *w);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", norm_2(res2->getDense() - prod(tmp3->getSym(), w->getSparse())) < tol, true);
  // Sparse * ...
  *res = prod(*tmp4, *vv);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", norm_2(res->getDense() - prod(tmp4->getSparse(), vv->getDense())) < tol, true);
  *res = prod(*tmp4, *ww);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", norm_2(res->getDense() - prod(tmp4->getSparse(), ww->getSparse())) < tol, true);
  // Triang * ...
  *res = prod(*tmp5, *v);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", norm_2(res->getDense() - prod(tmp5->getBanded(), v->getDense())) < tol, true);
  *res = prod(*tmp5, *w);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", norm_2(res->getDense() - prod(tmp5->getBanded(), w->getSparse())) < tol, true);
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  delete sv;
  delete v;
  delete w;
  delete vv;
  delete ww;
  delete res;
  delete res2;
  cout << "-->  test prod ended with success." << endl;
}

void SimpleMatrixTest::testProdBis()
{
  cout << "--> Test: prod. mat-vect (bis)" << endl;

  SiconosVector * y = new SimpleVector(size);
  SiconosVector * x = new SimpleVector(size, 4.3);
  SiconosVector * x1 = new SimpleVector(size - 2, 2.3);
  SiconosVector * x2 = new SimpleVector(2, 3.1);

  SiconosVector * xB = new BlockVector(x1, x2);
  SiconosVector * yB = new BlockVector(*xB);
  yB->zero();

  // Matrix - vector product

  // Simple = Simple * Simple
  prod(*A, *x, *y);
  double sum;
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*y)(i) - sum) < tol, true);
  }
  // Simple = Simple * Block
  prod(*A , *xB, *y);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*y)(i) - sum) < tol, true);
  }

  // Block = Simple * Simple
  prod(*A , *x, *yB);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*yB)(i) - sum) < tol, true);
  }

  // Block = Simple * Block
  prod(*A , *xB, *yB);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*yB)(i) - sum) < tol, true);
  }

  delete x;
  delete y;
  delete xB;
  delete yB;
  delete x1;
  delete x2;

  // Others or old stuff ...

  SiconosMatrix * tmp2 = new SimpleMatrix(*T);
  SiconosMatrix * tmp3 = new SimpleMatrix(*S);
  SiconosMatrix * tmp4 = new SimpleMatrix(*SP);
  SiconosMatrix * tmp5 = new SimpleMatrix(*Band2);
  SiconosVector * v = new SimpleVector(3);
  (*v)(0) = 1;
  (*v)(1) = 2;
  (*v)(2) = 3;
  SiconosVector * vv = new SimpleVector(4);
  (*vv)(0) = 1;
  (*vv)(1) = 2;
  (*vv)(2) = 3;
  SparseVect * sv = new SparseVect(3);
  (*sv)(0) = 4;
  (*sv)(1) = 5;
  (*sv)(2) = 6;
  SparseVect * sv2 = new SparseVect(4);
  (*sv2)(0) = 4;
  (*sv2)(1) = 5;
  (*sv2)(2) = 6;
  SiconosVector * w = new SimpleVector(*sv);
  SiconosVector * ww = new SimpleVector(*sv2);
  SiconosVector * res = new SimpleVector(4);
  SiconosVector * res2 = new SimpleVector(3);

  // Triang * ...
  prod(*tmp2, *v, *res2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", norm_2(res2->getDense() - prod(tmp2->getTriang(), v->getDense())) < tol, true);
  prod(*tmp2, *w, *res2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", norm_2(res2->getDense() - prod(tmp2->getTriang(), w->getSparse())) < tol, true);
  //   Sym * ...
  prod(*tmp3, *v, *res2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", norm_2(res2->getDense() - prod(tmp3->getSym(), v->getDense())) < tol, true);
  prod(*tmp3, *w, *res2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", norm_2(res2->getDense() - prod(tmp3->getSym(), w->getSparse())) < tol, true);
  // Sparse * ...
  prod(*tmp4, *vv, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", norm_2(res->getDense() - prod(tmp4->getSparse(), vv->getDense())) < tol, true);
  prod(*tmp4, *ww, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", norm_2(res->getDense() - prod(tmp4->getSparse(), ww->getSparse())) < tol, true);
  // Banded * ...
  prod(*tmp5, *v, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", norm_2(res->getDense() - prod(tmp5->getBanded(), v->getDense())) < tol, true);
  prod(*tmp5, *w, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", norm_2(res->getDense() - prod(tmp5->getBanded(), w->getSparse())) < tol, true);
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  delete sv;
  delete v;
  delete w;
  delete vv;
  delete ww;
  delete res;
  delete res2;
  cout << "-->  test prodBis ended with success." << endl;
}

void SimpleMatrixTest::testProdTer()
{
  cout << "--> Test: prod. mat-vect (ter)" << endl;

  SiconosVector * y = new SimpleVector(size);
  SiconosVector * x = new SimpleVector(size, 4.3);
  SiconosVector * x1 = new SimpleVector(size - 2, 2.3);
  SiconosVector * x2 = new SimpleVector(2, 3.1);

  SiconosVector * xB = new BlockVector(x1, x2);
  SiconosVector * yB = new BlockVector(*xB);
  yB->zero();

  // Matrix - vector product

  // Simple = Simple * Simple
  axpy_prod(*A, *x, *y, true);
  double sum;
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*y)(i) - sum) < tol, true);
  }

  SiconosVector * backUp = new SimpleVector(*y);
  // Simple += Simple * Simple
  axpy_prod(*A, *x, *y, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*y)(i) - sum - (*backUp)(i)) < tol, true);
  }

  // Simple = Simple * Block
  axpy_prod(*A , *xB, *y, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*y)(i) - sum) < tol, true);
  }

  *backUp = *y;
  // Simple += Simple * Block
  axpy_prod(*A , *xB, *y, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*y)(i) - sum - (*backUp)(i)) < tol, true);
  }

  // Block = Simple * Simple
  axpy_prod(*A , *x, *yB, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*yB)(i) - sum) < tol, true);
  }

  // Block += Simple * Simple
  *backUp = *yB;
  axpy_prod(*A , *x, *yB, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*yB)(i) - sum - (*backUp)(i)) < tol, true);
  }

  // Block = Simple * Block
  axpy_prod(*A , *xB, *yB, true);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*yB)(i) - sum) < tol, true);
  }

  // Block += Simple * Block
  *backUp = *yB;
  axpy_prod(*A , *xB, *yB, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*yB)(i) - sum - (*backUp)(i)) < tol, true);
  }

  delete backUp;
  delete x;
  delete y;
  delete xB;
  delete yB;
  delete x1;
  delete x2;

  // Others or old stuff ...

  SiconosMatrix * tmp2 = new SimpleMatrix(*T);
  SiconosMatrix * tmp3 = new SimpleMatrix(*S);
  SiconosMatrix * tmp4 = new SimpleMatrix(*SP);
  SiconosMatrix * tmp5 = new SimpleMatrix(*Band2);
  SiconosVector * v = new SimpleVector(3);
  (*v)(0) = 1;
  (*v)(1) = 2;
  (*v)(2) = 3;
  SiconosVector * vv = new SimpleVector(4);
  (*vv)(0) = 1;
  (*vv)(1) = 2;
  (*vv)(2) = 3;
  SparseVect * sv = new SparseVect(3);
  (*sv)(0) = 4;
  (*sv)(1) = 5;
  (*sv)(2) = 6;
  SparseVect * sv2 = new SparseVect(4);
  (*sv2)(0) = 4;
  (*sv2)(1) = 5;
  (*sv2)(2) = 6;
  SiconosVector * w = new SimpleVector(*sv);
  SiconosVector * ww = new SimpleVector(*sv2);
  SiconosVector * res = new SimpleVector(4);
  SiconosVector * res2 = new SimpleVector(3);

  // Triang * ...
  prod(*tmp2, *v, *res2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", norm_2(res2->getDense() - prod(tmp2->getTriang(), v->getDense())) < tol, true);
  prod(*tmp2, *w, *res2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", norm_2(res2->getDense() - prod(tmp2->getTriang(), w->getSparse())) < tol, true);
  //   Sym * ...
  prod(*tmp3, *v, *res2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", norm_2(res2->getDense() - prod(tmp3->getSym(), v->getDense())) < tol, true);
  prod(*tmp3, *w, *res2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", norm_2(res2->getDense() - prod(tmp3->getSym(), w->getSparse())) < tol, true);
  // Sparse * ...
  prod(*tmp4, *vv, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", norm_2(res->getDense() - prod(tmp4->getSparse(), vv->getDense())) < tol, true);
  prod(*tmp4, *ww, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", norm_2(res->getDense() - prod(tmp4->getSparse(), ww->getSparse())) < tol, true);
  // Banded * ...
  prod(*tmp5, *v, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", norm_2(res->getDense() - prod(tmp5->getBanded(), v->getDense())) < tol, true);
  prod(*tmp5, *w, *res);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", norm_2(res->getDense() - prod(tmp5->getBanded(), w->getSparse())) < tol, true);
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  delete sv;
  delete v;
  delete w;
  delete vv;
  delete ww;
  delete res;
  delete res2;
  cout << "-->  test prodTer ended with success." << endl;
}

void SimpleMatrixTest::testProd4() // y += A*x
{
  cout << "--> Test: prod. mat-vect (4)" << endl;

  SiconosVector * y = new SimpleVector(size);
  SiconosVector * x = new SimpleVector(size, 4.3);
  SiconosVector * x1 = new SimpleVector(size - 2, 2.3);
  SiconosVector * x2 = new SimpleVector(2, 3.1);

  SiconosVector * xB = new BlockVector(x1, x2);
  SiconosVector * yB = new BlockVector(*xB);
  yB->zero();

  // Matrix - vector product

  // Simple = Simple * Simple
  y->zero();
  prod(*A, *x, *y, false);
  prod(*A, *x, *y, false);
  double sum;
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += 2 * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*y)(i) - sum) < tol, true);
  }
  // Simple = Simple * Block
  y->zero();
  prod(*A , *xB, *y, false);
  prod(*A , *xB, *y, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += 2 * (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*y)(i) - sum) < tol, true);
  }

  // Block = Simple * Simple
  yB->zero();
  prod(*A , *x, *yB, false);
  prod(*A , *x, *yB, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += 2 * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*yB)(i) - sum) < tol, true);
  }

  // Block = Simple * Block
  yB->zero();
  prod(*A , *xB, *yB, false);
  prod(*A , *xB, *yB, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += 2 * (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*yB)(i) - sum) < tol, true);
  }

  delete x;
  delete y;
  delete xB;
  delete yB;
  delete x1;
  delete x2;
  cout << "-->  test prod4 ended with success." << endl;
}

void SimpleMatrixTest::testProd5() // y += a*A*x
{
  cout << "--> Test: prod. mat-vect (5)" << endl;

  SiconosVector * y = new SimpleVector(size);
  SiconosVector * x = new SimpleVector(size, 4.3);
  SiconosVector * x1 = new SimpleVector(size - 2, 2.3);
  SiconosVector * x2 = new SimpleVector(2, 3.1);

  SiconosVector * xB = new BlockVector(x1, x2);
  SiconosVector * yB = new BlockVector(*xB);
  yB->zero();

  // Matrix - vector product
  double a = 3.0;
  // Simple = Simple * Simple
  y->zero();
  prod(a, *A, *x, *y, false);
  prod(a, *A, *x, *y, false);
  double sum;
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += 2 * a * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*y)(i) - sum) < tol, true);
  }
  // Simple = Simple * Block
  y->zero();
  prod(a, *A , *xB, *y, false);
  prod(a, *A , *xB, *y, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += a * 2 * (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*y)(i) - sum) < tol, true);
  }

  // Block = Simple * Simple
  yB->zero();
  prod(a, *A , *x, *yB, false);
  prod(a, *A , *x, *yB, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += a * 2 * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*yB)(i) - sum) < tol, true);
  }

  // Block = Simple * Block
  yB->zero();
  prod(a, *A , *xB, *yB, false);
  prod(a, *A , *xB, *yB, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += a * 2 * (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*yB)(i) - sum) < tol, true);
  }

  delete x;
  delete y;
  delete xB;
  delete yB;
  delete x1;
  delete x2;
  cout << "-->  test prod5 ended with success." << endl;
}

void SimpleMatrixTest::testProd6() // y += trans(A)*x
{
  cout << "--> Test: prod. mat-vect (6)" << endl;

  SiconosVector * y = new SimpleVector(size);
  SiconosVector * x = new SimpleVector(size, 4.3);
  SiconosVector * x1 = new SimpleVector(size - 2, 2.3);
  SiconosVector * x2 = new SimpleVector(2, 3.1);

  SiconosVector * xB = new BlockVector(x1, x2);
  SiconosVector * yB = new BlockVector(*xB);
  yB->zero();

  SiconosMatrix * tmp = new SimpleMatrix(*A);
  tmp->trans();
  // Matrix - vector product

  // Simple = Simple * Simple
  y->zero();
  prod(*x, *A, *y);
  prod(*x, *A, *y, false);
  double sum;
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += 2 * (*tmp)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*y)(i) - sum) < tol, true);
  }
  // Simple = Simple * Block
  y->zero();
  prod(*xB, *A, *y);
  prod(*xB, *A, *y, false);

  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < size; ++j)
      sum += 2 * (*tmp)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*y)(i) - sum) < tol, true);
  }

  // Block = Simple * Simple
  yB->zero();
  prod(*x, *A , *yB);
  prod(*x, *A , *yB, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += 2 * (*tmp)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*yB)(i) - sum) < tol, true);
  }

  // Block = Simple * Block
  yB->zero();
  prod(*xB, *A , *yB);
  prod(*xB, *A , *yB, false);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += 2 * (*tmp)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*yB)(i) - sum) < tol, true);
  }

  delete tmp;
  delete x;
  delete y;
  delete xB;
  delete yB;
  delete x1;
  delete x2;
  cout << "-->  test prod6 ended with success." << endl;
}

void SimpleMatrixTest::testGemv()
{
  cout << "--> Test: gemv" << endl;

  SiconosVector * y = new SimpleVector(size, 1.0);
  SiconosVector * x = new SimpleVector(size, 4.3);

  SiconosVector * backUp = new SimpleVector(*y);

  gemv(*A, *x, *y);
  double sum;
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = 0.0;
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testgemv: ", fabs((*y)(i) - sum) < tol, true);
  }

  double a = 2.3;
  double b = 1.5;
  *y = *backUp;
  gemv(a, *A, *x, b, *y);

  for (unsigned int i = 0; i < size; ++i)
  {
    sum = b * (*backUp)(i);
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += a * (*A)(i, j) * (*x)(j) ;
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testgemv: ", fabs((*y)(i) - sum) < tol, true);
  }

  *y = *backUp;
  gemv(CblasNoTrans, a, *A, *x, b, *y);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = b * (*backUp)(i);
    for (unsigned int j = 0; j < A->size(1); ++j)
      sum += a * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testgemv: ", fabs((*y)(i) - sum) < tol, true);
  }
  *y = *backUp;
  gemv(CblasTrans, a, *A, *x, b, *y);
  for (unsigned int i = 0; i < size; ++i)
  {
    sum = b * (*backUp)(i);
    for (unsigned int j = 0; j < A->size(0); ++j)
      sum += a * (*A)(j, i) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testgemv: ", fabs((*y)(i) - sum) < tol, true);
  }


  delete x;
  delete y;
  delete backUp;

  cout << "-->  test gemv ended with success." << endl;
}

void SimpleMatrixTest::testGemm()
{
  cout << "--> Test: gemm." << endl;

  double a = 2.3;
  double b = 1.5;
  *C = *A;
  SiconosMatrix * backUp = new SimpleMatrix(*C);

  gemm(*A, *B, *C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGemm: ", norm_inf(*C->getDensePtr() - prod(*A->getDensePtr(), *B->getDensePtr())) < tol, true);

  *C = *backUp;
  gemm(a, *A, *B, b, *C);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGemm: ", norm_inf(*C->getDensePtr() - a * prod(*A->getDensePtr(), *B->getDensePtr()) - b**backUp->getDensePtr()) < tol, true);


  *C = *backUp;
  gemm(CblasNoTrans, CblasNoTrans, a, *A, *B, b, *C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGemm: ", norm_inf(*C->getDensePtr() - a * prod(*A->getDensePtr(), *B->getDensePtr()) - b**backUp->getDensePtr()) < tol, true);

  *C = *backUp;
  gemm(CblasTrans, CblasTrans, a, *A, *B, b, *C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGemm: ", norm_inf(*C->getDensePtr() - a * prod(trans(*A->getDensePtr()), trans(*B->getDensePtr())) - b**backUp->getDensePtr()) < tol, true);


  delete backUp;
  cout << "-->  test gemm ended with success." << endl;
}



void SimpleMatrixTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of SimpleMatrix Tests ===== " << endl;
  cout << "======================================" << endl;
}
