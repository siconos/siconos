/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
#ifndef __LagrangianDSXMLTest__
#define __LagrangianDSXMLTest__

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "LagrangianDSXML.h"
#include<vector>
#include<iostream>



class LagrangianDSXMLTest : public CppUnit::TestFixture
{


private:

  xmlDoc *doc;
  xmlNode *root;
  LagrangianDSXML ds;
  SiconosMatrix matrixRef;
  /*SiconosVector*/
  SimpleVector vectorRef;


  // on nomme la suite de tests
  CPPUNIT_TEST_SUITE(LagrangianDSXMLTest);

  // on ajoute les tests a effectuer

  // les tests qui doivent passer
  CPPUNIT_TEST(testGetNdof);
  CPPUNIT_TEST(testGetQVelocity);
  CPPUNIT_TEST(testGetMass);
  CPPUNIT_TEST(testIsPlugin);
  CPPUNIT_TEST(testGetPluginName);
  CPPUNIT_TEST(testGetMemory);

  //CPPUNIT_TEST_EXCEPTION(testIfTagIsNotPresent, XMLException);

  // on termine
  CPPUNIT_TEST_SUITE_END();

  // declaration de fonctions de test
  void testGetNdof();
  void testGetQVelocity();
  void testGetMemory();
  void testGetMass();
  void testIsPlugin();
  void testGetPluginName();

public:
  void setUp();
  void tearDown();
};

#endif
