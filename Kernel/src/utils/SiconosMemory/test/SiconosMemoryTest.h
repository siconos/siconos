//$Id: SiconosMemoryTest.h,v 1.4 2004/09/10 11:26:26 charlety Exp $
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
  vector<SiconosVector*> V;

public:
  void setUp();
  void tearDown();

};

#endif

//$Log: SiconosMemoryTest.h,v $
//Revision 1.4  2004/09/10 11:26:26  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.3  2004/07/27 14:41:16  charlety
//
//_ new tests for operator = of SiconosMemory
//_ modification of the makefile : the variable LIBS was false since the object
//  SiconosMemoryXML is linked with SiconosMemory.
//
//Revision 1.2  2004/07/08 09:13:43  charlety
//
//_ creation of a Siconos exception dedicated to the class SiconosMemory
//_ new tests for SiconosMemory
//
//Revision 1.1  2004/07/07 13:53:13  charlety
//
//_ First version of the memory object
//_ some little bugs corrected otherwhere
//