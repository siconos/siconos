//$id$

#ifndef SICONOSMEMORYXMLTEST_H
#define SICONOSMEMORYXMLTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "SiconosMemoryXML.h"

class SiconosMemoryXMLTest : public CppUnit::TestFixture
{
public:
  void setUp();
  void tearDown();

private:
  xmlDoc *doc;
  xmlNode *root;
  SiconosMemoryXML *smxml;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(SiconosMemoryXMLTest);

  // on ajoute les tests a effectuer

  // les tests qui doivent passer
  CPPUNIT_TEST(testHasMemory);

  // on termine
  CPPUNIT_TEST_SUITE_END();

  // declaration de fonctions de test
  void testHasMemory();
};

#endif // SICONOSMEMORYXMLTEST_H

//$log$