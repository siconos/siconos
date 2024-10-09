/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include <SiconosConfig.h>
#include <cppunit/extensions/HelperMacros.h>

#include "AVI.hpp"
#include "EulerMoreauOSI.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NormalConeNSL.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"

class AVITest : public CppUnit::TestFixture {
 private:
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

  unsigned int _n{2};
  double _h{0.1};
  double _t0{0.};
  double _T{10.};
  double _tol{1.e-12};
  double _theta{0.5};
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> _nsds{nullptr};
  std::shared_ptr<siconos::simulation::TimeStepping> _sim{nullptr};
  std::shared_ptr<siconos::modeling::FirstOrderLinearTIDS> _DS{nullptr};
  std::shared_ptr<siconos::simulation::TimeDiscretisation> _TD{nullptr};
  std::shared_ptr<siconos::integrators::OneStepIntegrator> _osi{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _A{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _b{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _x0{nullptr};

 public:
  AVITest() = default;
  void setUp();
  void tearDown();
};

#endif
