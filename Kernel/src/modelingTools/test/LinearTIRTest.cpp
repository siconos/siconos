#include "LinearTIRTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LinearTIRTest);


void LinearTIRTest::setUp()
{
  C = new SiconosMatrix("matC.dat", true);
  D = new SiconosMatrix("matD.dat", true);
  B = new SiconosMatrix("matB.dat", true);
  F = new SiconosMatrix("matF.dat", true);
  a = new SimpleVector(3);
  (*a)(0) = 4;
  (*a)(1) = 5;
  (*a)(2) = 6;
  e = new SimpleVector(2);
  (*e)(0) = 0.1;
  (*e)(1) = 0.1;
}

void LinearTIRTest::tearDown()
{
  delete e;
  delete a;
  delete F;
  delete B;
  delete C;
  delete D;
}

// default constructor

// Default  -> private -> no test
/*void LinearTIRTest::testBuildLinearTIR()
{
  LinearTIR * ltir = new LinearTIR();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSimpleVector : ", ltir->getType()==LINEARTIRELATION, true);
  delete ltir;
  cout << " Constructor LTIR 0 ok" << endl;
}
*/
// data constructor (1)
void LinearTIRTest::testBuildLinearTIR1()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIR : ", ltir->getB() == *B, true);
  delete ltir;
  cout << " Constructor LTIR 1 ok" << endl;

}

// data constructor (2)
void LinearTIRTest::testBuildLinearTIR2()
{
  LinearTIR * ltir = new LinearTIR(*C, *D, *F, *e, *B, *a);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRC : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRD : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRF : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRE : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRB : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLinearTIRA : ", ltir->getA() == *a, true);
  delete ltir;
  cout << " Constructor LTIR 2 ok" << endl;
}

// set C
void LinearTIRTest::testSetC()
{
  SiconosMatrix * tmp = new SiconosMatrix(*C);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*tmp, *B);
  ltir->setC(*C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", ltir->getC() == *C, true);
  delete ltir;
  delete tmp;
  cout << " testSetC ok" << endl;
}

// setCPtr
void LinearTIRTest::testSetCPtr()
{
  SiconosMatrix * tmp = new SiconosMatrix(*C);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*tmp, *B);
  ltir->setCPtr(C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", ltir->getC() == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", ltir->getCPtr() == C, true);
  delete ltir;
  delete tmp;
  cout << " test setCPtr ok" << endl;
}

// set D
void LinearTIRTest::testSetD()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setD(*D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetD: ", ltir->getD() == *D, true);
  delete ltir;
  cout << " test setD ok" << endl;
}

// setDPtr
void LinearTIRTest::testSetDPtr()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr : ", ltir->getD() == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", ltir->getDPtr() == D, true);
  delete ltir;
  cout << " test setDPtr ok" << endl;
}

// set F
void LinearTIRTest::testSetF()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setF(*F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetF: ", ltir->getF() == *F, true);
  delete ltir;
  cout << " test setF ok" << endl;
}

// setFPtr
void LinearTIRTest::testSetFPtr()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setFPtr(F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr : ", ltir->getF() == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", ltir->getFPtr() == F, true);
  delete ltir;
  cout << " test setFPtr ok" << endl;
}

// set E
void LinearTIRTest::testSetE()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setE(*e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetE: ", ltir->getE() == *e, true);
  delete ltir;
  cout << " test setE ok" << endl;
}

// setEPtr
void LinearTIRTest::testSetEPtr()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setEPtr(e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr : ", ltir->getE() == *e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", ltir->getEPtr() == e, true);
  delete ltir;
  cout << " test setEPtr ok" << endl;
}

// set B
void LinearTIRTest::testSetB()
{
  SiconosMatrix * tmp = new SiconosMatrix(*B);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*C, *tmp);
  ltir->setB(*B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetB: ", ltir->getB() == *B, true);
  delete ltir;
  delete tmp;
  cout << " test setB ok" << endl;
}

// setBPtr
void LinearTIRTest::testSetBPtr()
{
  SiconosMatrix * tmp = new SiconosMatrix(*B);
  tmp->zero();
  LinearTIR * ltir = new LinearTIR(*C, *tmp);
  ltir->setBPtr(B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", ltir->getB() == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", ltir->getBPtr() == B, true);
  delete ltir;
  delete tmp;
  cout << " test setBPtr ok" << endl;
}

// set A
void LinearTIRTest::testSetA()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setA(*a);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetA: ", ltir->getA() == *a, true);
  delete ltir;
  cout << " test setA ok" << endl;
}

// setAPtr
void LinearTIRTest::testSetAPtr()
{
  LinearTIR * ltir = new LinearTIR(*C, *B);
  ltir->setAPtr(a);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAPtr : ", ltir->getA() == *a, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAPtr: ", ltir->getAPtr() == a, true);
  delete ltir;
  cout << " test setAPtr ok" << endl;
}

void LinearTIRTest::End()
{
  cout << "====================================" << endl;
  cout << " ===== End of LinearTIR Tests ===== " << endl;
  cout << "====================================" << endl;
}
