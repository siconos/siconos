#ifndef __BoundaryConditionXMLTest__
#define __BoundaryConditionXMLTest__

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "LinearBCXML.h"
#include "PeriodicBCXML.h"
#include "NLinearBCXML.h"
#include "XMLException.h"



class BoundaryConditionXMLTest : public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root;
  LinearBCXML bcLinear;
  PeriodicBCXML bcPeriodic;
  NLinearBCXML bcNLinear;
  SiconosMatrix matrixRef;
  /*SiconosVector*/
  SimpleVector vectorRef;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(BoundaryConditionXMLTest);

  // on ajoute les tests a effectuer

  // les tests qui doivent passer

  CPPUNIT_TEST(testGetType);
  CPPUNIT_TEST(testGetOmega);

  //CPPUNIT_TEST_EXCEPTION(testIfTagIsNotPresent, XMLException);

  // on termine
  CPPUNIT_TEST_SUITE_END();

  // declaration de fonctions de test
  void testGetType();
  void testGetOmega();

public:
  void setUp();
  void tearDown();

};

#endif
