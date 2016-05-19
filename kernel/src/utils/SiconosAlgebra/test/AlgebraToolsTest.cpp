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
#include "SiconosConfig.h"
#include "AlgebraToolsTest.hpp"
#include "SimpleMatrix.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(AlgebraToolsTest);

//using namespace Siconos::algebra;


void AlgebraToolsTest::setUp()
{
  std::string fic1 = "matForExp.dat"; // 3X3
  
  A.reset(new SimpleMatrix(fic1, true));
}

void AlgebraToolsTest::tearDown()
{}

//______________________________________________________________________________

void AlgebraToolsTest::testExpm() 
{
  std::cout << "====================================" <<std::endl;
  std::cout << "=== AlgebraTools tests start ...=== " <<std::endl;
  std::cout << "====================================" <<std::endl;
  std::cout << "--> Test: expm " <<std::endl;
  std::string fic2 = "matSolForExp.dat"; // 3X3
  SP::SimpleMatrix ref(new SimpleMatrix(fic2, true));
  SP::SimpleMatrix Exp(new SimpleMatrix(3,3));
  
  Siconos::algebra::tools::expm(*A, *Exp);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testExpm : ", (*ref - *Exp).normInf() < 1e-6, true);
  std::cout << "--> Expm test ended with success." <<std::endl;
}


void AlgebraToolsTest::End()
{
  std::cout << "======================================" <<std::endl;
  std::cout << " ===== End of AlgebraTools Tests ===== " <<std::endl;
  std::cout << "======================================" <<std::endl;
}
