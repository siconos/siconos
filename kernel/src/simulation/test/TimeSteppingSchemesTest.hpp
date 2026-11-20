/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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
#ifndef TimeSteppingSchemesTest_hpp
#define TimeSteppingSchemesTest_hpp

#include <cppunit/extensions/HelperMacros.h>

/** Integration (non-regression) tests of the time-stepping schemes.

    Each test builds a small nonsmooth mechanical (or electrical) system, integrates it over
    a short time with one scheme, and compares 10 samples of the trajectory (time, state,
    gaps and multipliers) with a reference file `<test name>.ref` of this directory.

    The systems cover: linear and nonlinear Lagrangian systems, Newton-Euler bodies, single
    and multiple contacts, simultaneous impacts, resting contacts, Coulomb friction (sliding
    and sticking), bilateral constraints, a free rigid body and a first-order system.

    To create or update the reference files after an intended change of the results, run the
    test with the environment variable SICONOS_REGENERATE_TEST_REFS set: the references are
    then written into the source directory (and the test passes).
*/
class TimeSteppingSchemesTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(TimeSteppingSchemesTest);

  // Column of beads (linear Lagrangian, 4 unilateral contacts, multiple impacts)
  CPPUNIT_TEST(beadsMoreauJean);
  CPPUNIT_TEST(beadsCombinedProjection);
  CPPUNIT_TEST(beadsDirectProjection);
  CPPUNIT_TEST(beadsD1MinusAccelerationLevel);
  CPPUNIT_TEST(beadsD1MinusVelocityLevel);
  CPPUNIT_TEST(beadsSchatzmanPaoli);

  // Linear oscillators against obstacles (diagonal Lagrangian system)
  CPPUNIT_TEST(oscillatorsMoreauJeanBilbao);

  // Disk with Coulomb friction in 2D (sliding then rolling)
  CPPUNIT_TEST(disk2DFrictionMoreauJean);

  // Sphere with Coulomb friction in 3D (local and global formulations)
  CPPUNIT_TEST(sphere3DFrictionMoreauJean);
  CPPUNIT_TEST(sphere3DFrictionMoreauJeanGOSI);

  // Rocking block (nonlinear relations, 2 contacts)
  CPPUNIT_TEST(rockingBlockMoreauJean);
  CPPUNIT_TEST(rockingBlockCombinedProjection);
  CPPUNIT_TEST(rockingBlockCombinedProjectionFromOrigin);
  CPPUNIT_TEST(rockingBlockDirectProjection);
  CPPUNIT_TEST(rockingBlockD1Minus);

  // Double pendulum (nonlinear Lagrangian system, nonlinear relations, 2 contacts)
  CPPUNIT_TEST(doublePendulumMoreauJean);
  CPPUNIT_TEST(doublePendulumD1Minus);

  // Pendulum in Cartesian coordinates (bilateral constraint + unilateral contact)
  CPPUNIT_TEST(bilateralPendulumMoreauJean);
  CPPUNIT_TEST(bilateralPendulumCombinedProjection);

  // Newton-Euler cube landing on a plane (8 corners)
  CPPUNIT_TEST(cubeFrictionMoreauJean);
  CPPUNIT_TEST(cubeFrictionMoreauJeanGOSI);
  CPPUNIT_TEST(cubeMoreauJean);
  CPPUNIT_TEST(cubeDirectProjection);
  CPPUNIT_TEST(cubeCombinedProjection);

  // Free Newton-Euler body (gyroscopic effects, intermediate axis)
  CPPUNIT_TEST(freeRigidBodyMoreauJean);

  // First-order system: RLC circuit with a diode
  CPPUNIT_TEST(circuitRLCDEulerMoreau);
  CPPUNIT_TEST(circuitRLCDEulerMoreauResiduals);

  // Variants added for the coverage of the time-stepping code
  CPPUNIT_TEST(beadsManagedMoreauJean);
  CPPUNIT_TEST(beadsMoreauJeanSparseBlockStorage);
  CPPUNIT_TEST(sphere3DFrictionMoreauJeanReducedDirect);
  CPPUNIT_TEST(springChainClampedMoreauJean);
  CPPUNIT_TEST(springChainMoreauJean);
  CPPUNIT_TEST(frictionOscillatorMoreauJean);
  CPPUNIT_TEST(frictionOscillatorMoreauJeanGOSI);
  CPPUNIT_TEST(frictionOscillatorSparseMoreauJean);
  CPPUNIT_TEST(frictionOscillatorSparseMoreauJeanGOSI);
  CPPUNIT_TEST(sphere3DFrictionMoreauJeanGOSIGlobalReduced);
  CPPUNIT_TEST(bilateralPendulumDirectProjection);
  CPPUNIT_TEST(doublePendulumMoreauJeanNonlinearFull);
  CPPUNIT_TEST(stackMoreauJean);
  CPPUNIT_TEST(stackDirectProjection);
  CPPUNIT_TEST(stackMoreauJeanSymmetricGraph);
  CPPUNIT_TEST(springChainSchatzmanPaoli);

  CPPUNIT_TEST_SUITE_END();

 public:
  void beadsMoreauJean();
  void beadsCombinedProjection();
  void beadsDirectProjection();
  void beadsD1MinusAccelerationLevel();
  void beadsD1MinusVelocityLevel();
  void beadsSchatzmanPaoli();
  void oscillatorsMoreauJeanBilbao();
  void disk2DFrictionMoreauJean();
  void sphere3DFrictionMoreauJean();
  void sphere3DFrictionMoreauJeanGOSI();
  void rockingBlockMoreauJean();
  void rockingBlockCombinedProjection();
  void rockingBlockCombinedProjectionFromOrigin();
  void rockingBlockDirectProjection();
  void rockingBlockD1Minus();
  void doublePendulumMoreauJean();
  void doublePendulumD1Minus();
  void bilateralPendulumMoreauJean();
  void bilateralPendulumCombinedProjection();
  void cubeFrictionMoreauJean();
  void cubeFrictionMoreauJeanGOSI();
  void cubeMoreauJean();
  void cubeDirectProjection();
  void cubeCombinedProjection();
  void freeRigidBodyMoreauJean();
  void circuitRLCDEulerMoreau();
  void circuitRLCDEulerMoreauResiduals();

 private:
  void circuitRLCD(bool residuals);

 public:
  void beadsManagedMoreauJean();
  void beadsMoreauJeanSparseBlockStorage();
  void sphere3DFrictionMoreauJeanReducedDirect();
  void springChainClampedMoreauJean();
  void springChainMoreauJean();
  void frictionOscillatorMoreauJean();
  void frictionOscillatorMoreauJeanGOSI();
  void frictionOscillatorSparseMoreauJean();
  void frictionOscillatorSparseMoreauJeanGOSI();
  void sphere3DFrictionMoreauJeanGOSIGlobalReduced();
  void bilateralPendulumDirectProjection();
  void doublePendulumMoreauJeanNonlinearFull();
  void stackMoreauJean();
  void stackDirectProjection();
  void stackMoreauJeanSymmetricGraph();
  void springChainSchatzmanPaoli();
};

#endif
