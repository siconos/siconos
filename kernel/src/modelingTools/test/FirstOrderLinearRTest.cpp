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
#include "FirstOrderLinearRTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearRTest);


void FirstOrderLinearRTest::setUp()
{
  C.reset(new SimpleMatrix("matC.dat", true));
  D.reset(new SimpleMatrix("matD.dat", true));
  B.reset(new SimpleMatrix("matB.dat", true));
  F.reset(new SimpleMatrix("matF.dat", true));
  e.reset(new SiconosVector(1));
  //   Cp.reset(new FirstOrderLinearR::PluggedMatrix("TestPlugin:C"));
  //   Dp.reset(new FirstOrderLinearR::PluggedMatrix("TestPlugin:D"));
  //   Bp.reset(new FirstOrderLinearR::PluggedMatrix("TestPlugin:B"));
  //   Fp.reset(new FirstOrderLinearR::PluggedMatrix("TestPlugin:F"));
  //  ep.reset(new Plugged_Vector_FTime("TestPlugin:e"));
  (*e)(0) = 0.1;
}

void FirstOrderLinearRTest::tearDown()
{}

// data constructor (1)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR1()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::FirstOrderLinearR folr(new FirstOrderLinearR("TestPlugin:C", "TestPlugin:B"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->getSubType() == RELATION::LinearR, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->C()->pluginName()=="TestPlugin:C", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR1 : ", folr->B()->pluginName()=="TestPlugin:B", true);
  std::cout << "--> Constructor 1 test ended with success." <<std::endl;
}

void FirstOrderLinearRTest::testBuildFirstOrderLinearR3()
{
  std::cout << "--> Test: constructor 3." <<std::endl;

  SP::FirstOrderLinearR folr(new FirstOrderLinearR("TestPlugin:C", "TestPlugin:D", "TestPlugin:F", "TestPlugin:e", "TestPlugin:B"));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->getSubType() == RELATION::LinearR, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->C()->pluginName()=="TestPlugin:C", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->D()->pluginName()=="TestPlugin:D", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->F()->pluginName()=="TestPlugin:F", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->e()->pluginName()=="TestPlugin:e", true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR3 : ", folr->B()->pluginName()=="TestPlugin:B", true);
  std::cout << "--> Constructor 3 test ended with success." <<std::endl;
}

// data constructor (4)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR4()
{
  std::cout << "--> Test: constructor 4." <<std::endl;
  SP::FirstOrderLinearR folr(new FirstOrderLinearR(C, B));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4b : ", folr->B() == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4c : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR4d : ", folr->getSubType() == RELATION::LinearR, true);
  std::cout << "--> Constructor 4 test ended with success." <<std::endl;
}

// data constructor (5)
void FirstOrderLinearRTest::testBuildFirstOrderLinearR5()
{
  std::cout << "--> Test: constructor 5." <<std::endl;
  SP::FirstOrderLinearR folr(new FirstOrderLinearR(C, D, F, e, B));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5b : ", folr->D() == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5c : ", folr->F() == F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5d : ", folr->e() == e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5e : ", folr->B() == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5f : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearR5g : ", folr->getSubType() == RELATION::LinearR, true);
  std::cout << "--> Constructor 5 test ended with success." <<std::endl;
}

// set C as a matrix and then plug it

// setCPtr
void FirstOrderLinearRTest::testSetCPtr()
{
  //   std::cout << "--> Test: setCPtr." <<std::endl;
  //   FirstOrderLinearR::SP_PluggedMatrix tmp(new FirstOrderLinearR::PluggedMatrix(*C));
  //   tmp->zero();
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*tmp,*B));
  //   folr->setCPtr(Cp);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->getC()==*Cp, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->C()==Cp, true);
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->C()->isPlugged(), true);
  // //   folr->setComputeCFunction("TestPlugin.so","C");
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->C()->isPlugged()==true, true);
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC : ", folr->C()->pluginName()=="TestPlugin:C", true);

  std::cout << "--> setCPtr test ended with success." <<std::endl;
}
// set C as a plugin
void FirstOrderLinearRTest::testSetCPtr2()
{
  //   std::cout << "--> Test: setCPtr2." <<std::endl;
  //   FirstOrderLinearR::SP_PluggedMatrix tmp(new FirstOrderLinearR::PluggedMatrix("TestPlugin:D"));
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(tmp,Bp));
  //   folr->setCPtr(Cp);
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC2 : ", folr->C()->isPlugged()==true, true);
  // //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetC2 : ", folr->C()->pluginName()=="TestPlugin:C", true);

  //   std::cout << "--> setCPtr2 test ended with success." <<std::endl;
}

// set D

// setDPtr
void FirstOrderLinearRTest::testSetDPtr()
{
  //   std::cout << "--> Test: setDPtr." <<std::endl;
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*B));
  //   folr->setDPtr(Dp);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr : ", folr->getD()==*Dp, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", folr->D()==Dp, true);
  //   std::cout << "--> setDPtr test ended with success." <<std::endl;
}

// set F

// setFPtr
void FirstOrderLinearRTest::testSetFPtr()
{
  //   std::cout << "--> Test: setFPtr." <<std::endl;
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*B));
  //   folr->setFPtr(Fp);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr : ", folr->getF()==*Fp, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", folr->F()==Fp, true);
  //   std::cout << "--> setFPtr test ended with success." <<std::endl;
}

// set E


// setEPtr
void FirstOrderLinearRTest::testSetEPtr()
{
  //   std::cout << "--> Test: setEPtr." <<std::endl;
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*B));
  //   folr->setEPtr(ep);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr : ", folr->getE()==*ep, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", folr->e()==ep, true);
  //   std::cout << "--> setEPtr test ended with success." <<std::endl;
}

// set B

// setBPtr
void FirstOrderLinearRTest::testSetBPtr()
{
  //   std::cout << "--> Test: setBPtr." <<std::endl;
  //   FirstOrderLinearR::SP_PluggedMatrix tmp(new FirstOrderLinearR::PluggedMatrix(*B));
  //   tmp->zero();
  //   SP::FirstOrderLinearR folr(new FirstOrderLinearR(*C,*tmp));
  //   folr->setBPtr(Bp);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", folr->getB()==*Bp, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", folr->B()==Bp, true);
  //   std::cout << "--> setBPtr test ended with success." <<std::endl;
}
