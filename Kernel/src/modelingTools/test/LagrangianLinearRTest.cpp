#include "LagrangianLinearRTest.h"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianLinearRTest);


void LagrangianLinearRTest::setUp()
{
  H = new SiconosMatrix("matH.dat", true);
  b = new SimpleVector(2);
  (*b)(0) = 4;
  (*b)(1) = 5;
}

void LagrangianLinearRTest::tearDown()
{
  delete b;
  delete H;
}

// data constructor (1)
void LagrangianLinearRTest::testBuildLagrangianLinearR1()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H, *b);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR : ", llr->getH() == *H, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR : ", llr->getB() == *b, true);
  delete llr;
  cout << " Constructor LLR 1 ok" << endl;
}

// data constructor (2)
void LagrangianLinearRTest::testBuildLagrangianLinearR2()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearRC : ", llr->getH() == *H, true);
  delete llr;
  cout << " Constructor LLR 2 ok" << endl;
}

// copy constructor
void LagrangianLinearRTest::testBuildLagrangianLinearR3()
{
  Relation * rel = new LagrangianLinearR(*H, *b);
  LagrangianLinearR *llr = new LagrangianLinearR(*rel);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR : ", llr->getH() == *H, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearR : ", llr->getB() == *b, true);
  delete rel;
  delete llr;
  cout << " Constructor LLR 3 ok" << endl;
}

// set H
void LagrangianLinearRTest::testSetH()
{
  SiconosMatrix * tmp = new SiconosMatrix(*H);
  tmp->zero();
  LagrangianLinearR * llr = new LagrangianLinearR(*tmp, *b);
  llr->setH(*H);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetH : ", llr->getH() == *H, true);
  delete llr;
  delete tmp;
  cout << " testSetH ok" << endl;
}

// setHPtr
void LagrangianLinearRTest::testSetHPtr()
{
  SiconosMatrix * tmp = new SiconosMatrix(*H);
  tmp->zero();
  LagrangianLinearR * llr = new LagrangianLinearR(*tmp, *b);
  llr->setHPtr(H);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetHPtr : ", llr->getH() == *H, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetHPtr : ", llr->getHPtr() == H, true);
  delete llr;
  delete tmp;
  cout << " test setHPtr ok" << endl;
}

// setB
void LagrangianLinearRTest::testSetB()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H);
  llr->setB(*b);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetb: ", llr->getB() == *b, true);
  delete llr;
  cout << " test setB ok" << endl;
}

// setBPtr
void LagrangianLinearRTest::testSetBPtr()
{
  LagrangianLinearR * llr = new LagrangianLinearR(*H);
  llr->setBPtr(b);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", llr->getB() == *b, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", llr->getBPtr() == b, true);
  delete llr;
  cout << " test setBPtr ok" << endl;
}
