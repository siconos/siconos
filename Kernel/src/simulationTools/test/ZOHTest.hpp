/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef __ZOHTest__
#define __ZOHTest__

#include <cppunit/extensions/HelperMacros.h>
#include "FirstOrderLinearTIDS.hpp"
#include "ZeroOrderHold.hpp"
#include "Model.hpp"
#include "TimeStepping.hpp"
#include "TimeDiscretisation.hpp"
#include "ioMatrix.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "RelayNSL.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Relay.hpp"

class ZOHTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(ZOHTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(ZOHTest);

  // tests to be done ...

  CPPUNIT_TEST(testMatrixExp0);
  CPPUNIT_TEST(testMatrixExp1);
  CPPUNIT_TEST(testMatrixIntegration1);
  CPPUNIT_TEST(testMatrixIntegration2);
  CPPUNIT_TEST(testMatrixIntegration3);
  CPPUNIT_TEST(testMatrixIntegration4);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test
  void init();
  void testMatrixExp0();
  void testMatrixExp1();
  void testMatrixIntegration1();
  void testMatrixIntegration2();
  void testMatrixIntegration3();
  void testMatrixIntegration4();
  // Members

  unsigned int _n;
  double _h;
  double _t0;
  double _T;
  double _tol;
  SP::Model _model;
  SP::TimeStepping _sim;
  SP::FirstOrderLinearDS _DS;
  SP::TimeDiscretisation _TD;
  SP::ZeroOrderHold _ZOH;
  SP::SiconosMatrix _A;
  SP::SiconosVector _b;
  SP::SiconosVector _x0;


public:
  void setUp();
  void tearDown();

};

#endif




