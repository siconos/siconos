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
#include "FirstOrderLinearRTest.hpp"

#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearRTest);

void FirstOrderLinearRTest::setUp()
{
  C = std::make_shared<siconos::algebra::SiconosMatrix>(siconos::algebra::readMatrixFromFile("matC.dat"));
  D = std::make_shared<siconos::algebra::SiconosMatrix>(siconos::algebra::readMatrixFromFile("matD.dat"));
  B = std::make_shared<siconos::algebra::SiconosMatrix>(siconos::algebra::readMatrixFromFile("matB.dat"));
  F = std::make_shared<siconos::algebra::SiconosMatrix>(siconos::algebra::readMatrixFromFile("matF.dat"));



  e = std::make_shared<siconos::algebra::SiconosVector>(1);
  (*e)(0) = 0.1;
}

void FirstOrderLinearRTest::tearDown() {}

// data constructor (1)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR1()
{
  std::cout << "--> Test: constructor 1." << std::endl;
  auto folr =
      std::make_shared<siconos::modeling::FirstOrderLinearR>("TestPlugin:C", "TestPlugin:B");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ",
                               folr->getType() == siconos::modeling::RelationType::FirstOrder,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderLinearR1 : ",
      folr->getSubType() == siconos::modeling::RelationSubType::LinearR, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ",
  //   folr->C()->pluginName()=="TestPlugin:C", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ",
  //   folr->B()->pluginName()=="TestPlugin:B", true);
  std::cout << "--> Constructor 1 test ended with success." << std::endl;
}

void FirstOrderLinearRTest::testBuildFirstOrderLinearR3()
{
  std::cout << "--> Test: constructor 3." << std::endl;

  auto folr = std::make_shared<siconos::modeling::FirstOrderLinearR>(
      "TestPlugin:C", "TestPlugin:D", "TestPlugin:F", "TestPlugin:e", "TestPlugin:B");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ",
                               folr->getType() == siconos::modeling::RelationType::FirstOrder,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderLinearR3 : ",
      folr->getSubType() == siconos::modeling::RelationSubType::LinearR, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ",
  //   folr->C()->pluginName()=="TestPlugin:C", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ",
  //   folr->D()->pluginName()=="TestPlugin:D", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ",
  //   folr->F()->pluginName()=="TestPlugin:F", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ",
  //   folr->e()->pluginName()=="TestPlugin:e", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ",
  //   folr->B()->pluginName()=="TestPlugin:B", true);
  std::cout << "--> Constructor 3 test ended with success." << std::endl;
}

// data constructor (4)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR4()
{
  std::cout << "--> Test: constructor 4." << std::endl;
  auto folr = std::make_shared<siconos::modeling::FirstOrderLinearR>(C, B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4b : ", folr->B() == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4c : ",
                               folr->getType() == siconos::modeling::RelationType::FirstOrder,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderLinearR4d : ",
      folr->getSubType() == siconos::modeling::RelationSubType::LinearR, true);
  std::cout << "--> Constructor 4 test ended with success." << std::endl;
}

// data constructor (5)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR5()
{
  std::cout << "--> Test: constructor 5." << std::endl;
  auto folr = std::make_shared<siconos::modeling::FirstOrderLinearR>(C, D, F, e, B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5b : ", folr->D() == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5c : ", folr->F() == F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5d : ", folr->e() == e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5e : ", folr->B() == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5f : ",
                               folr->getType() == siconos::modeling::RelationType::FirstOrder,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildFirstOrderLinearR5g : ",
      folr->getSubType() == siconos::modeling::RelationSubType::LinearR, true);
  std::cout << "--> Constructor 5 test ended with success." << std::endl;
}

// // set C as a matrix and then plug it

// // setCPtr
// void FirstOrderLinearRTest::testSetCPtr()
// {
//   //   std::cout << "--> Test: setCPtr." <<std::endl;
//   //   FirstOrderLinearR::SP_PluggedMatrix tmp=
//   std::make_shared<siconos::modeling::FirstOrderLinearR>::PluggedMatrix(*C));
//   //   tmp->zero();
//   //   auto folr= std::make_shared<siconos::modeling::FirstOrderLinearR>(*tmp,*B);
//   //   folr->setCPtr(Cp);
//   //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->getC()==*Cp, true);
//   //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->C()==Cp, true);
//   // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->C()->isPlugged(), true);
//   // //   folr->setComputeCFunction("TestPlugin.so","C");
//   // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->C()->isPlugged()==true, true);
//   // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ",
//   folr->C()->pluginName()=="TestPlugin:C", true);

//   std::cout << "--> setCPtr test ended with success." <<std::endl;
// }
// // set C as a plugin
// void FirstOrderLinearRTest::testSetCPtr2()
// {
//   //   std::cout << "--> Test: setCPtr2." <<std::endl;
//   //   FirstOrderLinearR::SP_PluggedMatrix tmp=
//   std::make_shared<siconos::modeling::FirstOrderLinearR>::PluggedMatrix("TestPlugin:D");
//   //   auto folr= std::make_shared<siconos::modeling::FirstOrderLinearR>(tmp,Bp));
//   //   folr->setCPtr(Cp);
//   // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC2 : ", folr->C()->isPlugged()==true, true);
//   // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC2 : ",
//   folr->C()->pluginName()=="TestPlugin:C", true);

//   //   std::cout << "--> setCPtr2 test ended with success." <<std::endl;
// }

// // set D

// // setDPtr
// void FirstOrderLinearRTest::testSetDPtr()
// {
//   //   std::cout << "--> Test: setDPtr." <<std::endl;
//   //   auto folr= std::make_shared<siconos::modeling::FirstOrderLinearR>(*C,*B);
//   //   folr->setDPtr(Dp);
//   //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr : ", folr->getD()==*Dp, true);
//   //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", folr->D()==Dp, true);
//   //   std::cout << "--> setDPtr test ended with success." <<std::endl;
// }

// // set F

// // setFPtr
// void FirstOrderLinearRTest::testSetFPtr()
// {
//   //   std::cout << "--> Test: setFPtr." <<std::endl;
//   //   auto folr= std::make_shared<siconos::modeling::FirstOrderLinearR>(*C,*B);
//   //   folr->setFPtr(Fp);
//   //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr : ", folr->getF()==*Fp, true);
//   //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", folr->F()==Fp, true);
//   //   std::cout << "--> setFPtr test ended with success." <<std::endl;
// }

// // set E

// // setEPtr
// void FirstOrderLinearRTest::testSetEPtr()
// {
//   //   std::cout << "--> Test: setEPtr." <<std::endl;
//   //   auto folr= std::make_shared<siconos::modeling::FirstOrderLinearR>(*C,*B);
//   //   folr->setEPtr(ep);
//   //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr : ", folr->getE()==*ep, true);
//   //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", folr->e()==ep, true);
//   //   std::cout << "--> setEPtr test ended with success." <<std::endl;
// }

// // set B

// // setBPtr
// void FirstOrderLinearRTest::testSetBPtr()
// {
//   //   std::cout << "--> Test: setBPtr." <<std::endl;
//   //   FirstOrderLinearR::SP_PluggedMatrix tmp=
//   std::make_shared<siconos::modeling::FirstOrderLinearR>::PluggedMatrix(*B);
//   //   tmp->zero();
//   //   auto folr= std::make_shared<siconos::modeling::FirstOrderLinearR>(*C,*tmp));
//   //   folr->setBPtr(Bp);
//   //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", folr->getB()==*Bp, true);
//   //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", folr->B()==Bp, true);
//   //   std::cout << "--> setBPtr test ended with success." <<std::endl;
// }
