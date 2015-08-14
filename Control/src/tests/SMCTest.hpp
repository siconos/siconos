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
#ifndef __SMCTest__
#define __SMCTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SiconosFwd.hpp"
#include "SiconosControlFwd.hpp"
#include <FirstOrderLinearTIDS.hpp>

#include <NumericsConfig.h>

class SMCTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SMCTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SMCTest);

  // tests to be done ...

  CPPUNIT_TEST(test_iSMC_ZOH);
  CPPUNIT_TEST(test_iSMC_Lsodar);
  CPPUNIT_TEST(test_eSMC_ZOH);
  CPPUNIT_TEST(test_eSMC_Lsodar);
#ifdef HAS_EXTREME_POINT_ALGO
  CPPUNIT_TEST(test_itw_ZOH);
  CPPUNIT_TEST(test_itw_Lsodar);
#endif

  CPPUNIT_TEST_SUITE_END();

  void init();
  void init2();
  void initTwisting();
  void test_iSMC_ZOH();
  void test_iSMC_Lsodar();
  void test_eSMC_ZOH();
  void test_eSMC_Lsodar();

#ifdef HAS_EXTREME_POINT_ALGO
  void test_itw_ZOH();
  void test_itw_Lsodar();
#endif
  // Members

  unsigned int _n;
  double _h;
  double _t0;
  double _T;
  double _tol;
  double _beta;
  double _xFinal;
  SP::FirstOrderLinearTIDS _DS;
  SP::SiconosMatrix _A;
  SP::SimpleMatrix _B;
  SP::SimpleMatrix _C;
  SP::SimpleMatrix _Csurface;
  SP::SiconosVector _b;
  SP::SiconosVector _x0;
  SP::SiconosVector _K;
  SP::LinearSensor _sensor;
  SP::LinearSMC _iSMC;
  SP::ExplicitLinearSMC _eSMC;
#ifdef HAS_EXTREME_POINT_ALGO
  SP::Twisting _itw;
#endif


public:

  SMCTest(): _n(2), _h(0.05), _t0(0.0), _T(100.0), _tol(7e-11), _beta(0.1) {}
  void setUp();
  void tearDown();

};

#endif




