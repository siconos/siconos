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
#include "SolidLinearTIDSTest.hpp"

#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "MeshUtils.hpp"
#include "FiniteElementModel.hpp"
#include "Material.hpp"
#include "SolidLinearTIDS.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SolidLinearTIDSTest);

void SolidLinearTIDSTest::setUp() {

  mesh = siconos::mechanics::fem::createMeshFromGMSH2("square_6.msh");
  int bulk_material_tag = 1;
  double density = 7800.;
  auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(density, 210e9, 1 / 3.);
  materials = {{bulk_material_tag, mat1}};

  std::shared_ptr<siconos::algebra::SiconosMatrix> denseS;
  denseS = std::make_shared<siconos::algebra::SimpleMatrix>("S.dat", true);
  S = std::make_shared<siconos::algebra::SimpleMatrix>(12,12,siconos::algebra::UblasType::SPARSE);
  for (int i = 0; i< 4; i++){
    for (int j = 0; j< 3; j++){
      for (int k = 0; k< 3; k++){
        S->setValue(3*i+j,3*i+k,denseS->getValue(3*i+j,3*i+k));
      }
    }
  }

  std::shared_ptr<siconos::algebra::SiconosMatrix> denseB;
  denseB = std::make_shared<siconos::algebra::SimpleMatrix>("B.dat", true);
  B = std::make_shared<siconos::algebra::SimpleMatrix>(12,10,siconos::algebra::UblasType::SPARSE);
  for (int i = 0; i< 12; i++){
    for (int j = 0; j< 10; j++){
      B->setValue(i,j,denseB->getValue(i,j));
    }
  }


}

void SolidLinearTIDSTest::tearDown() {}

// Mass, K, C
void SolidLinearTIDSTest::testBuildSolidLinearTIDS1() {
  std::cout << "--> Test: constructor 1." << std::endl;
  auto FEsolid = std::make_shared<siconos::mechanics::fem::SolidLinearTIDS>(
      mesh, materials, siconos::algebra::UblasType::SPARSE);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSolidLinearTIDS1 : ", FEsolid->dimension() == 10,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSolidLinearTIDS1 : ", FEsolid->velocityDimension() == 10,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSolidLinearTIDS1 : ", FEsolid->stressDimension() == 12,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSolidLinearTIDS1 : ", FEsolid->FEModel()->elements().size() == 4,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSolidLinearTIDS1 : ", (*FEsolid->stress())(0) == 0.0,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testS : ",((*FEsolid->S())-*S).normInf() < 1.0e-16,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testB: ",((*FEsolid->B())-*B).normInf() < 1.0e-16,
                               true);
  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}
