/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#ifndef __AVITest__
#define __AVITest__

#include <cppunit/extensions/HelperMacros.h>
#include "FirstOrderLinearTIDS.hpp"
#include "TimeStepping.hpp"
#include "TimeDiscretisation.hpp"
#include "ioMatrix.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "NormalConeNSL.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "AVI.hpp"
#include "EulerMoreauOSI.hpp"

#include <SiconosConfig.h>

class AVITest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(AVITest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(AVITest);

  // tests to be done ...

#ifdef HAS_EXTREME_POINT_ALGO
  CPPUNIT_TEST(testAVI);
#endif

  CPPUNIT_TEST_SUITE_END();

  void init();
  void testAVI();

  unsigned int _n;
  double _h;
  double _t0;
  double _T;
  double _tol;
  double _theta;
  SP::NonSmoothDynamicalSystem _nsds;
  SP::TimeStepping _sim;
  SP::FirstOrderLinearTIDS _DS;
  SP::TimeDiscretisation _TD;
  SP::OneStepIntegrator _osi;
  SP::SiconosMatrix _A;
  SP::SiconosVector _b;
  SP::SiconosVector _x0;


public:

  AVITest(): _n(2), _h(0.1), _t0(0.0), _T(10.0), _tol(1e-12), _theta(0.5) {}
  void setUp();
  void tearDown();

};

#endif




