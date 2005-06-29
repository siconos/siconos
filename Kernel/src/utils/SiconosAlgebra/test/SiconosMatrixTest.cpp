#include "SiconosMatrixTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMatrixTest);


void SiconosMatrixTest::setUp()
{
  int i, j;
  int Arow = 10;
  int Acol = 10;
  int Brow = 10;
  int Bcol = 5;
  int Vsize = 10;

  A = SiconosMatrix(Arow, Acol);
  B = SiconosMatrix(Brow, Bcol);
  /*C==A*/
  C = SiconosMatrix(Arow, Acol);


  /* init A, C */
  srand((unsigned)time(NULL));
  for (i = 0; i < Arow; i++)
    for (int j = 0; j < Acol; j++)
    {
      A(i, j) = rand() % 10 + 20;
      C(i, j) = A(i, j);
    }

  /* init B */
  for (i = 0; i < Brow; i++)
    for (j = 0; j < Bcol; j++)
      B(i, j) = rand() % 100 - 50;


  /* init SV */
  vector<double> vtmp(Vsize);
  for (i = 0; i < Vsize; i++)
    vtmp.at(i) = rand() % 10 + 20;

  SV.setValues(vtmp);
}

void SiconosMatrixTest::tearDown()
{ }

//______________________________________________________________________________

void SiconosMatrixTest::testConstructor1()
{
  unsigned int row = 120;
  unsigned int col = 500;
  SiconosMatrix X(row, col);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : X.size(0) == row", X.size(0), row);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : X.size(1) == col", X.size(1), col);
  cout << " SiconosMatrixTest >>> testConstructor1 .................................................. OK\n ";
}

void SiconosMatrixTest::testConstructor2()
{

  LaVectorDouble LVD(10);
  for (int i = 0; i < LVD.size(); i++)
    LVD(i) = i;

  SiconosMatrix X(LVD, 5, 2);

  SiconosMatrix Xtrue(5, 2);
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 2; j++)
      Xtrue(i, j) = 2 * i + j;

  CPPUNIT_ASSERT_MESSAGE("testConstructor2 : X == Xtrue", X == Xtrue);
  cout << "SiconosMatrixTest >>> testConstructor2 .................................................. OK\n ";
}

void SiconosMatrixTest::testEquality()
{
  CPPUNIT_ASSERT_MESSAGE("testEquality : A==A", A == A);
  CPPUNIT_ASSERT_MESSAGE("testEquality : B==B", B == B);
  CPPUNIT_ASSERT_MESSAGE("testEquality : C==C", C == C);

  CPPUNIT_ASSERT_MESSAGE("testEquality : A!=B", A != B);
  CPPUNIT_ASSERT_MESSAGE("testEquality : B!=C", B != C);
  CPPUNIT_ASSERT_MESSAGE("testEquality : A==C", A == C);
  C(5, 5) = C(5, 5) * 2;
  CPPUNIT_ASSERT_MESSAGE("testEquality : A!=C", A != C);
  cout << "SiconosMatrixTest >>> testEquality ...................................................... OK\n ";
}


void SiconosMatrixTest::testLinearSolve()
{
  SiconosMatrix X;
  SiconosMatrix Xbefore = X;
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : Xbefore==X", Xbefore == X);
  X = A.linearSolve(B);
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : A==A", A == A);
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : B==B", B == B);
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : X!=A", X != A);
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : X!=B", X != B);
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : Xbefore!=X", Xbefore != X);
  cout << "SiconosMatrixTest >>> testLinearSolve ................................................... OK\n ";
}

void SiconosMatrixTest::testReadWriteBinary()
{
  SiconosMatrix X = B;
  A.write("totoMat.bin", "binary");
  X.read("totoMat.bin", "binary");

  CPPUNIT_ASSERT_MESSAGE("testReadWriteBinary : A == X ", A == X);
  SiconosMatrix Y("totoMat.bin", false);
  CPPUNIT_ASSERT_MESSAGE("testReadWriteBinary : A == Y ", A == Y);
  cout << "SiconosMatrixTest >>> testReadWriteBinary ............................................... OK\n ";
}

void SiconosMatrixTest::testReadWriteAscii()
{
  SiconosMatrix X = B;
  A.write("totoMat.ascii", "ascii");
  X.read("totoMat.ascii", "ascii");

  CPPUNIT_ASSERT_MESSAGE("testReadWriteAscii : A == X ", A == X);
  SiconosMatrix Y("totoMat.ascii", true);
  CPPUNIT_ASSERT_MESSAGE("testReadWriteAscii : A == Y ", A == Y);
  cout << "SiconosMatrixTest >>> testReadWriteAscii ................................................ OK\n ";
}

void SiconosMatrixTest::testAffectation()
{
  SiconosMatrix X = A;
  SiconosMatrix Y;
  Y = A;
  CPPUNIT_ASSERT_MESSAGE("testAffectation : A == X ", A == X);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : &A!=&X ", &A != &X);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : A == Y ", A == Y);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : &A!=&Y ", &A != &Y);
  A(1, 1) = A(1, 1) * 2;
  CPPUNIT_ASSERT_MESSAGE("testAffectation : A != X ", A != X);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : A != Y ", A != Y);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : X == Y ", X == Y);

  cout << "SiconosMatrixTest >>> testAffectation ................................................... OK\n ";
}

void SiconosMatrixTest::testOperator()
{
  cout << endl;
  B = A + A;
  B = A - A;
  B = A * 2;
  B = 2 * A;
  B = A / 2;
  B = A ^ 0;
  B = A ^ 5;
  B = A.multTranspose(A);
  B = A.multTranspose(B);

  CPPUNIT_ASSERT_MESSAGE("testOperator : A+A == A*2 ", A + A == A * 2);
  CPPUNIT_ASSERT_MESSAGE("testOperator : A+A == 2*A ", A + A == 2 * A);
  CPPUNIT_ASSERT_MESSAGE("testOperator : A^1 == A ", (A ^ 1) == A);
  CPPUNIT_ASSERT_MESSAGE("testOperator : A^2 == A*A ", (A ^ 2) == A * A);
  CPPUNIT_ASSERT_MESSAGE("testOperator : A^3 == A*A*A ", (A ^ 3) == A * A * A);
  CPPUNIT_ASSERT_MESSAGE("testOperator : A^10 == A*A*A*A*A*A*A*A*A*A ", (A ^ 10) == A * A * A * A * A * A * A * A * A * A);

  cout << "SiconosMatrixTest >>> testOperator ...................................................... OK\n ";
}

void SiconosMatrixTest::testAddRow()
{
  SiconosMatrix X = A;
  /*SiconosVector*/
  SimpleVector V2 = X.getRow(2);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : X == A ", X == A);
  X.addRow(2, SV);
  /*SiconosVector*/
  SimpleVector V3 = X.getRow(2);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : X != A ", X != A);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : V2 != V3 ", V2 != V3);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : V3 == SV ", V3 == SV);

  cout << "SiconosMatrixTest >>> testAddRow ........................................................ OK\n ";
}

void SiconosMatrixTest::testSizeException()
{
  A + B;
}

void SiconosMatrixTest::testConstructorException()
{
  LaVectorDouble LVD(10);
  for (int i = 0; i < LVD.size(); i++)
    LVD(i) = i;
  SiconosMatrix X(LVD, 4, 2);
}

void SiconosMatrixTest::testBlockMatrixCopy1()
{
  SiconosMatrix A(1, 2), B(3, 3), C(3, 3);
  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  B(0, 0) = 0.0;
  B(0, 1) = 0.0;
  B(0, 2) = 0.0;
  B(1, 0) = 0.0;
  B(1, 1) = 0.0;
  B(1, 2) = 0.0;
  B(2, 0) = 0.0;
  B(2, 1) = 0.0;
  B(2, 2) = 0.0;
  C(0, 0) = 0.0;
  C(0, 1) = 0.0;
  C(0, 2) = 0.0;
  C(1, 0) = 1.0;
  C(1, 1) = 2.0;
  C(1, 2) = 0.0;
  C(2, 0) = 0.0;
  C(2, 1) = 0.0;
  C(2, 2) = 0.0;

  B.blockMatrixCopy(A, 1, 0);
  CPPUNIT_ASSERT_MESSAGE("testBlockMatrixCopy1 : B.blockMatrixCopy( A, x, y) ", B == C);
  cout << "SiconosMatrixTest >>> testBlockMatrixCopy1 ........................................................ OK\n ";
}

void SiconosMatrixTest::testBlockMatrixCopy2()
{
  SiconosMatrix A(1, 1), B(1, 1), C(1, 1);
  A(0, 0) = 1.0;
  B(0, 0) = 0.0;
  C(0, 0) = 1.0;

  B.blockMatrixCopy(A, 0, 0);
  CPPUNIT_ASSERT_MESSAGE("testBlockMatrixCopy2 : B.blockMatrixCopy( A, x, y) ", B == C);
  cout << "SiconosMatrixTest >>> testBlockMatrixCopy2 ........................................................ OK\n ";
}

void SiconosMatrixTest::testBlockMatrixCopy3()
{
  SiconosMatrix A(2, 2), B(4, 4), C(4, 4);
  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  A(1, 0) = 3.0;
  A(1, 1) = 4.0;
  B(0, 0) = 0.0;
  B(0, 1) = 0.0;
  B(0, 2) = 0.0;
  B(0, 3) = 0.0;
  B(1, 0) = 0.0;
  B(1, 1) = 0.0;
  B(1, 2) = 0.0;
  B(1, 3) = 0.0;
  B(2, 0) = 0.0;
  B(2, 1) = 0.0;
  B(2, 2) = 0.0;
  B(2, 3) = 0.0;
  B(3, 0) = 0.0;
  B(3, 1) = 0.0;
  B(3, 2) = 0.0;
  B(3, 3) = 0.0;
  C(0, 0) = 1.0;
  C(0, 1) = 2.0;
  C(0, 2) = 0.0;
  C(0, 3) = 0.0;
  C(1, 0) = 3.0;
  C(1, 1) = 4.0;
  C(1, 2) = 0.0;
  C(1, 3) = 0.0;
  C(2, 0) = 0.0;
  C(2, 1) = 0.0;
  C(2, 2) = 0.0;
  C(2, 3) = 0.0;
  C(3, 0) = 0.0;
  C(3, 1) = 0.0;
  C(3, 2) = 0.0;
  C(3, 3) = 0.0;

  B.blockMatrixCopy(A, 0, 0);
  CPPUNIT_ASSERT_MESSAGE("testBlockMatrixCopy3 : B.blockMatrixCopy( A, x, y) ", B == C);
  cout << "SiconosMatrixTest >>> testBlockMatrixCopy3 ........................................................ OK\n ";
}

void SiconosMatrixTest::testBlockMatrixCopy4()
{
  SiconosMatrix A(2, 2), B(4, 4), C(4, 4);
  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  A(1, 0) = 3.0;
  A(1, 1) = 4.0;
  B(0, 0) = 0.0;
  B(0, 1) = 0.0;
  B(0, 2) = 0.0;
  B(0, 3) = 0.0;
  B(1, 0) = 0.0;
  B(1, 1) = 0.0;
  B(1, 2) = 0.0;
  B(1, 3) = 0.0;
  B(2, 0) = 0.0;
  B(2, 1) = 0.0;
  B(2, 2) = 0.0;
  B(2, 3) = 0.0;
  B(3, 0) = 0.0;
  B(3, 1) = 0.0;
  B(3, 2) = 0.0;
  B(3, 3) = 0.0;
  C(0, 0) = 0.0;
  C(0, 1) = 0.0;
  C(0, 2) = 0.0;
  C(0, 3) = 0.0;
  C(1, 0) = 0.0;
  C(1, 1) = 0.0;
  C(1, 2) = 0.0;
  C(1, 3) = 0.0;
  C(2, 0) = 0.0;
  C(2, 1) = 0.0;
  C(2, 2) = 1.0;
  C(2, 3) = 2.0;
  C(3, 0) = 0.0;
  C(3, 1) = 0.0;
  C(3, 2) = 3.0;
  C(3, 3) = 4.0;

  B.blockMatrixCopy(A, 2, 2);
  CPPUNIT_ASSERT_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", B == C);
  cout << "SiconosMatrixTest >>> testBlockMatrixCopy4 ........................................................ OK\n ";
}

void SiconosMatrixTest::testBlockMatrixCopyException1()
{
  SiconosMatrix A(2, 2), B(4, 4), C(4, 4);
  B.blockMatrixCopy(A, 3, 3);
}

void SiconosMatrixTest::testBlockMatrixCopyException2()
{
  SiconosMatrix A(6, 6), B(4, 4), C(4, 4);
  B.blockMatrixCopy(A, 0, 0);
}
