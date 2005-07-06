#ifndef __SiconosMemoryTest__
#define __SiconosMemoryTest__

#include <cppunit/extensions/HelperMacros.h>
#include "CompositeVector.h"
#include "SiconosMemory.h"
#include "SiconosMemoryException.h"
#include <deque>


#define SIZE 10
#define M_SIZE 3

using namespace std;

class SiconosMemoryTest : public CppUnit::TestFixture
{


private:

  // Test suite
  CPPUNIT_TEST_SUITE(SiconosMemoryTest);

  CPPUNIT_TEST(testBuildMemory1);
  CPPUNIT_TEST(testBuildMemory2);
  CPPUNIT_TEST(testBuildMemory3);
  CPPUNIT_TEST(testBuildMemory4);
  CPPUNIT_TEST(testSetVectorMemory);
  CPPUNIT_TEST(testGetSiconosVector);
  CPPUNIT_TEST(testSwap);
  CPPUNIT_TEST(testOperatorEqual);

  CPPUNIT_TEST_EXCEPTION(testMemoryException, SiconosMemoryException);
  CPPUNIT_TEST_EXCEPTION(testMemoryException1, SiconosMemoryException);

  //CPPUNIT_TEST_FAIL(testFail);

  CPPUNIT_TEST_SUITE_END();

  void testBuildMemory1();
  void testBuildMemory2();
  void testBuildMemory3();
  void testBuildMemory4();
  void testSetVectorMemory();
  void testGetSiconosVector();
  void testSwap();
  void testOperatorEqual();

  void testMemoryException();
  void testMemoryException1();
  //void testFail();

  deque<SiconosVector*> V1;
  deque<SiconosVector*> V2;
  deque<SiconosVector*> V3;
  SiconosVector * q1;
  SiconosVector * q2;
  SiconosVector * q3;
  SiconosVector *c1;
  SiconosVector *c2;
  unsigned int sizeMem;
public:
  void setUp();
  void tearDown();

};

#endif

