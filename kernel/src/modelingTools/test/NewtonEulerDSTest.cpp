/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "NewtonEulerDSTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(NewtonEulerDSTest);


void NewtonEulerDSTest::setUp()
{
  q0.reset(new SiconosVector(7));
  (*q0)(0) = 1;
  (*q0)(1) = 2;
  (*q0)(2) = 3;
  (*q0)(3) = 0.0;
  (*q0)(4) = 1.0;
  (*q0)(5) = 0.0;
  (*q0)(6) = 0.0;
  velocity0.reset(new SiconosVector(6));
  (*velocity0)(0) = 4;
  (*velocity0)(1) = 5;
  (*velocity0)(2) = 6;
  (*velocity0)(3) = 7;
  (*velocity0)(4) = 8;
  (*velocity0)(5) = 9;

  inertia.reset(new SimpleMatrix(3, 3));
  (*inertia)(0, 0) = 1;
  (*inertia)(1, 1) = 2;
  (*inertia)(2, 2) = 3;

  mass = 10.0;
}


void NewtonEulerDSTest::tearDown()
{}

// constructor from data
void NewtonEulerDSTest::testBuildNewtonEulerDS1()
{
  std::cout << "--> Test: constructor 1." <<std::endl;

  SP::NewtonEulerDS ds(new NewtonEulerDS(q0, velocity0, mass,  inertia ));
  double time = 1.5;
  ds->initialize(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1A : ", Type::value(*ds) == Type::NewtonEulerDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1B : ", ds->number() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", ds->dimension() == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", ds->getqDim() == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", ds->scalarMass() == mass, true);

  SP::SimpleMatrix massMatrix(new SimpleMatrix(6,6));
  massMatrix->setValue(0, 0, mass);
  massMatrix->setValue(1, 1, mass);
  massMatrix->setValue(2, 2, mass);

  Index dimIndex(2);
  dimIndex[0] = 3;
  dimIndex[1] = 3;
  Index startIndex(4);
  startIndex[0] = 0;
  startIndex[1] = 0;
  startIndex[2] = 3;
  startIndex[3] = 3;
  setBlock(inertia, massMatrix, dimIndex, startIndex);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", *(ds->mass()) == *(massMatrix), true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNewtonEulerDS1D : ", ds->computeKineticEnergy() == 595.0, true);

  std::cout << "--> Constructor 1 test ended with success." <<std::endl;
}

// void NewtonEulerDSTest::testcomputeDS()
// {
//   std::cout << "-->Test: computeDS." <<std::endl;
//   DynamicalSystem * ds(new NewtonEulerDS(tmpxml2));
//   SP::NewtonEulerDS copy =  std11::static_pointer_cast<NewtonEulerDS>(ds);
//   double time = 1.5;
//   ds->initialize("EventDriven", time);
//   ds->computeRhs(time);
//   std::cout << "-->Test: computeDS." <<std::endl;
//   ds->computeJacobianRhsx(time);
//   std::cout << "-->Test: computeDS." <<std::endl;
//   SimpleMatrix M(3, 3);
//   M(0, 0) = 1;
//   M(1, 1) = 2;
//   M(2, 2) = 3;
//   SP::SiconosMatrix jx = ds->jacobianRhsx();
//   SP::SiconosVector vf = ds->rhs();

//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSI : ", *(vf->vector(0)) == *velocity0, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSJ : ", prod(M, *(vf->vector(1))) == (copy->getFExt() - copy->getFInt() - copy->getFGyr()) , true);

//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(M, *(jx->block(1, 0))) == (copy->getJacobianFL(0)) , true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(M, *(jx->block(1, 1))) == (copy->getJacobianFL(1)) , true);
//   std::cout << "--> computeDS test ended with success." <<std::endl;

// }

void NewtonEulerDSTest::End()
{
  std::cout << "======================================" <<std::endl;
  std::cout << " ===== End of NewtonEulerDS tests ===== " <<std::endl;
  std::cout << "======================================" <<std::endl;
}
