/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#ifndef __EulerMoreauTest__
#define __EulerMoreauTest__

#include <cppunit/extensions/HelperMacros.h>
#include "FirstOrderLinearTIDS.hpp"
#include "EulerMoreauOSI.hpp"
#include "TimeStepping.hpp"
#include "TimeDiscretisation.hpp"
#include "ioMatrix.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "RelayNSL.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Relay.hpp"

class EulerMoreauTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(EulerMoreauTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(EulerMoreauTest);

  // tests to be done ...

  CPPUNIT_TEST(testCstGradTIDS);
  CPPUNIT_TEST(testCstGradDS);
  CPPUNIT_TEST(testCstGradNLDS);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test
  void init(bool initDS);
  void testCstGradTIDS();
  void testCstGradDS();
  void testCstGradNLDS();
  // Members

  unsigned int _n;
  double _h;
  double _t0;
  double _T;
  double _tol;
  SP::NonSmoothDynamicalSystem _model;
  SP::TimeStepping _sim;
  SP::DynamicalSystem _DS;
  SP::TimeDiscretisation _TD;
  SP::EulerMoreauOSI _EulerMoreau;
  SP::SiconosMatrix _A;
  SP::SiconosVector _b;
  SP::SiconosVector _x0;


public:

  EulerMoreauTest(): _n(2), _h(0.1), _t0(0.0), _T(10.0), _tol(1e-12) {}
  void setUp();
  void tearDown();

};

#endif




