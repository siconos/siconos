#ifndef NUMERICSTEST_H
#define NUMERICSTEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "SiconosNumerics.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
#include <fstream>


using namespace std;

class TestNumerics : public CppUnit::TestFixture
{
private:
  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(TestNumerics);

  // on ajoute les tests a effectuer :
  // les tests qui doivent passer
  CPPUNIT_TEST(testLcp);
  //  CPPUNIT_TEST(testLcp2);
  CPPUNIT_TEST(testRp);
  CPPUNIT_TEST(testCfp);
  CPPUNIT_TEST(testCfd);

  CPPUNIT_TEST(testDLSODE);
  CPPUNIT_TEST(testDLSODES);
  CPPUNIT_TEST(testDLSODA);
  CPPUNIT_TEST(testDLSODAR);
  CPPUNIT_TEST(testDLSODPK);
  CPPUNIT_TEST(testDLSODKR);
  CPPUNIT_TEST(testDLSODI);
  CPPUNIT_TEST(testDLSOIBT);
  CPPUNIT_TEST(testDLSODIS);

  // on termine
  CPPUNIT_TEST_SUITE_END();

  //  // declaration de fonctions de test
  void testLcp();
  void testLcp2();
  void testRp();
  void testCfd();
  void testCfp();

  void testDLSODE();
  void testDLSODES();
  void testDLSODA();
  void testDLSODAR();
  void testDLSODPK();
  void testDLSODKR();
  void testDLSODI();
  void testDLSOIBT();
  void testDLSODIS();

  // declararion des variables de tests


public:
  TestNumerics();
  ~TestNumerics();

  void setUp();
  void tearDown();
};

#endif
