/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#ifndef __DynamicalSystemXMLTest__
#define __DynamicalSystemXMLTest__

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "DynamicalSystemXML.h"
#include<iostream>
#include<vector>


class DynamicalSystemXMLTest : public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root;
  DynamicalSystemXML ds;
  SiconosMatrix matrixRef;
  /*SiconosVector*/
  SimpleVector vectorRef;

  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(DynamicalSystemXMLTest);

  // on ajoute les tests a effectuer

  // les tests qui doivent passer
  CPPUNIT_TEST(testGetNumber);
  CPPUNIT_TEST(testGetType);
  CPPUNIT_TEST(testGetId);
  CPPUNIT_TEST(testGetN);
  CPPUNIT_TEST(testGetStepsInMemory);
  CPPUNIT_TEST(testGetX);
  CPPUNIT_TEST(testGetPluginName);

  //CPPUNIT_TEST_EXCEPTION(testIfTagIsNotPresent, XMLException);

  // on termine
  CPPUNIT_TEST_SUITE_END();

  // declaration de fonctions de test
  void testGetNumber();
  void testGetType();
  void testGetId();
  void testGetN();
  void testGetStepsInMemory();
  void testGetX();
  void testGetPluginName();

public:
  void setUp();
  void tearDown();

};

#endif
