/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
//$id$

#ifndef SICONOSMEMORYXMLTEST_H
#define SICONOSMEMORYXMLTEST_H

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <cstdlib>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "SiconosMemoryXML.hpp"
#include <iostream>

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

