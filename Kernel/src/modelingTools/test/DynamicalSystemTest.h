#ifndef __DynamicalSystemTest__
#define __DynamicalSystemTest__

#include <cppunit/extensions/HelperMacros.h>
#include "DynamicalSystem.h"
#include "RuntimeException.h"
#include "XMLException.h"

class DynamicalSystemTest : public CppUnit::TestFixture
{

private:

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(DynamicalSystemTest);

  // tests to be done ...

  CPPUNIT_TEST(testBuildDynamicalSystem1);
  CPPUNIT_TEST(testBuildDynamicalSystem2);
  CPPUNIT_TEST(testBuildDynamicalSystem3);
  CPPUNIT_TEST(testBuildDynamicalSystem4);
  CPPUNIT_TEST(testSetX0);
  CPPUNIT_TEST(testSetX0Ptr);
  CPPUNIT_TEST_EXCEPTION(testSetX02, RuntimeException);
  CPPUNIT_TEST(testSetX);
  CPPUNIT_TEST(testSetXPtr);
  CPPUNIT_TEST_EXCEPTION(testSetX2, RuntimeException);
  CPPUNIT_TEST(testSetXDot);
  CPPUNIT_TEST(testSetXDotPtr);
  CPPUNIT_TEST_EXCEPTION(testSetXDot2, RuntimeException);
  CPPUNIT_TEST(testSetXFree);
  CPPUNIT_TEST(testSetXFreePtr);
  CPPUNIT_TEST_EXCEPTION(testSetXFree2, RuntimeException);
  CPPUNIT_TEST(testSetR);
  CPPUNIT_TEST(testSetRPtr);
  CPPUNIT_TEST_EXCEPTION(testSetR2, RuntimeException);
  CPPUNIT_TEST(testSetJacobianX);
  CPPUNIT_TEST(testSetJacobianXPtr);
  CPPUNIT_TEST_EXCEPTION(testSetJacobianX2, RuntimeException);
  CPPUNIT_TEST(testSetU);
  CPPUNIT_TEST(testSetUPtr);
  CPPUNIT_TEST_EXCEPTION(testSetU2, RuntimeException);
  CPPUNIT_TEST(testSetT);
  CPPUNIT_TEST(testSetTPtr);
  CPPUNIT_TEST_EXCEPTION(testSetT2, RuntimeException);
  CPPUNIT_TEST(testInitMemory);
  CPPUNIT_TEST(testSwap);
  CPPUNIT_TEST(End);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test

  void testBuildDynamicalSystem1();
  void testBuildDynamicalSystem2();
  void testBuildDynamicalSystem3();
  void testBuildDynamicalSystem4();
  void testSetX0();
  void testSetX0Ptr();
  void testSetX02();
  void testSetX();
  void testSetXPtr();
  void testSetX2();
  void testSetXDot();
  void testSetXDotPtr();
  void testSetXDot2();
  void testSetXFree();
  void testSetXFreePtr();
  void testSetXFree2();
  void testSetR();
  void testSetRPtr();
  void testSetR2();
  void testSetJacobianX();
  void testSetJacobianXPtr();
  void testSetJacobianX2();
  void testSetU();
  void testSetUPtr();
  void testSetU2();
  void testSetT();
  void testSetTPtr();
  void testSetT2();
  void testInitMemory();
  void testSwap();
  void End();

  // Members

  SimpleVector * x0, *u0;
  SiconosMatrix *T0, *J0;
  xmlNode * node1 , *node2;
  DynamicalSystemXML* tmpxml1, * tmpxml2;
public:
  void setUp();
  void tearDown();

};

#endif




