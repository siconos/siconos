#ifndef __SiconosMemoryTest__
#define __SiconosMemoryTest__

#include <cppunit/extensions/HelperMacros.h>
#include <vector>
#include "NewSiconosVector.h"
#include "SiconosMemory.h"
#include "SiconosMemoryException.h"


#define SIZE 10
#define M_SIZE 3

using namespace std;

class SiconosMemoryTest : public CppUnit::TestFixture
{


private:

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(SiconosMemoryTest);

  // on ajoute les tests a effectuer :

  // les tests qui doivent passer
  CPPUNIT_TEST(testBuildMemory);
  CPPUNIT_TEST(testBuildMemory1);
  CPPUNIT_TEST(testBuildMemory2);
  CPPUNIT_TEST(testBuildMemory3);
  CPPUNIT_TEST(testGetVectorMemory);
  CPPUNIT_TEST(testGetSiconosVector);
  CPPUNIT_TEST(testGetSiconosVector1);
  CPPUNIT_TEST(testSetVectorMemory);
  CPPUNIT_TEST(testSetVectorMemory1);
  CPPUNIT_TEST(testSetVectorMemory2);
  CPPUNIT_TEST(testSwap);
  CPPUNIT_TEST(testOperatorEqual);
  CPPUNIT_TEST(testOperatorEqual1);

  CPPUNIT_TEST_EXCEPTION(testMemoryException, SiconosMemoryException);
  CPPUNIT_TEST_EXCEPTION(testMemoryException1, SiconosMemoryException);
  CPPUNIT_TEST_EXCEPTION(testMemoryException2, SiconosMemoryException);
  CPPUNIT_TEST_EXCEPTION(testMemoryException3, SiconosMemoryException);

  // les tests qui doivent echouer
  //CPPUNIT_TEST_FAIL(testFail);

  // on termine
  CPPUNIT_TEST_SUITE_END();

  // declaration de fonctions de test
  void testBuildMemory();
  void testBuildMemory1();
  void testBuildMemory2();
  void testBuildMemory3();
  void testGetVectorMemory();
  void testGetSiconosVector();
  void testGetSiconosVector1();
  void testSetVectorMemory();
  void testSetVectorMemory1();
  void testSetVectorMemory2();
  void testSwap();
  void testOperatorEqual();
  void testOperatorEqual1();

  void testMemoryException();
  void testMemoryException1();
  void testMemoryException2();
  void testMemoryException3();
  void testFail();

  // declaration des variables de tests

  //SiconosVector *A, *B, *C;
  vector<SiconosVector*> V1;
  vector<SiconosVector*> V2;
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

