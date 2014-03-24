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
#include "LagrangianDSTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianDSTest);


void LagrangianDSTest::setUp()
{
  q0.reset(new SiconosVector(3));
  (*q0)(0) = 1;
  (*q0)(1) = 2;
  (*q0)(2) = 3;
  velocity0.reset(new SiconosVector(3));
  (*velocity0)(0) = 4;
  (*velocity0)(1) = 5;
  (*velocity0)(2) = 6;

  u0.reset(new SiconosVector(2));
  (*u0)(0) = 4;
  (*u0)(1) = 5;

  mass.reset(new SimpleMatrix(3, 3));
  (*mass)(0, 0) = 1;
  (*mass)(1, 1) = 2;
  (*mass)(2, 2) = 3;



void LagrangianDSTest::tearDown()
{}

// constructor from data
void LagrangianDSTest::testBuildLagrangianDS4()
{
  std::cout << "--> Test: constructor 4." <<std::endl;

  SP::LagrangianDS ds(new LagrangianDS(13, *q0, *velocity0, (*mass)));
  double time = 1.5;
  ds->initialize("TimeStepping", time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4A : ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4B : ", ds->number() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4D : ", ds->getNdof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4L : ", ds->getMass() == (*mass), true);

  map<string, bool> isPl = ds->getIsPlugin();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4Q : ", isPl["fExt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4R : ", isPl["fInt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4S : ", isPl["NNL"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4T : ", isPl["jacobianFIntq"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4T : ", isPl["jacobianFIntVelocity"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4T : ", isPl["jacobianNNLq"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4T : ", isPl["jacobianNNLVelocity"], false);
  std::cout << "--> Constructor 4 test ended with success." <<std::endl;
}

// constructor from data
void LagrangianDSTest::testBuildLagrangianDS5()
{
  std::cout << "--> Test: constructor 5." <<std::endl;
  std::string plugin = "TestPlugin:computeMass";
  SP::DynamicalSystem ds(new LagrangianDS(13, *q0, *velocity0, plugin));
  double time = 1.5;
  ds->initialize("TimeStepping", time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5A : ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5B : ", ds->number() == 13, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5D : ", std11::static_pointer_cast<LagrangianDS>(ds)->getNdof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5L : ", std11::static_pointer_cast<LagrangianDS>(ds)->getMass() == (*mass), true);

  map<string, bool> isPl = ds->getIsPlugin();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5M : ", isPl["mass"], true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5Q : ", isPl["fExt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5R : ", isPl["fInt"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5S : ", isPl["NNL"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5T : ", isPl["jacobianFIntq"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5T : ", isPl["jacobianFIntVelocity"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5T : ", isPl["jacobianNNLq"], false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5T : ", isPl["jacobianNNLVelocity"], false);
  std::cout << "--> Constructor 5 test ended with success." <<std::endl;
}

void LagrangianDSTest::testcomputeDS()
{
  std::cout << "-->Test: computeDS." <<std::endl;
  DynamicalSystem * ds(new LagrangianDS(tmpxml2));
  SP::LagrangianDS copy =  std11::static_pointer_cast<LagrangianDS>(ds);
  double time = 1.5;
  ds->initialize("EventDriven", time);
  ds->computeRhs(time);
  std::cout << "-->Test: computeDS." <<std::endl;
  ds->computeJacobianRhsx(time);
  std::cout << "-->Test: computeDS." <<std::endl;
  SimpleMatrix M(3, 3);
  M(0, 0) = 1;
  M(1, 1) = 2;
  M(2, 2) = 3;
  SP::SiconosMatrix jx = ds->jacobianRhsx();
  SP::SiconosVector vf = ds->rhs();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSI : ", *(vf->vector(0)) == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSJ : ", prod(M, *(vf->vector(1))) == (copy->getFExt() - copy->getFInt() - copy->getNNL()) , true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(M, *(jx->block(1, 0))) == (copy->getJacobianFL(0)) , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(M, *(jx->block(1, 1))) == (copy->getJacobianFL(1)) , true);
  std::cout << "--> computeDS test ended with success." <<std::endl;


}
void LagrangianDSTest::End()
{
  std::cout << "======================================" <<std::endl;
  std::cout << " ===== End of LagrangianDS tests ===== " <<std::endl;
  std::cout << "======================================" <<std::endl;
}
