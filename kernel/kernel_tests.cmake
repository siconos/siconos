include(tools4tests)

if(WITH_${COMPONENT}_TESTING)
  
  find_package(CPPUNIT REQUIRED)

  # Main test driver for cppunit tests
  set(TEST_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/test/TestMain.cpp)
  set(SIMPLE_TEST_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/SiconosMemory/test/TestMain.cpp)

  # ---- Siconos Algebra tests ----
  begin_tests(src/utils/SiconosAlgebra/test)

  new_test_1(
    NAME testSiconosAlgebra
    SOURCES BlockMatrixTest.cpp  SimpleMatrixTest.cpp BlockVectorTest.cpp  SiconosVectorTest.cpp EigenProblemsTest.cpp AlgebraToolsTest.cpp TestMain.cpp
    DEPS "numerics;CPPUNIT::CPPUNIT;externals"
    )

  # ---- Siconos Memory tests ----
  begin_tests(src/utils/SiconosMemory/test)

  new_test_1(
    NAME testSiconosMemory
    SOURCES SiconosMemoryTest.cpp TestMain.cpp
    DEPS "numerics;CPPUNIT::CPPUNIT"
    )

  # ---- Siconos tools tests ----
  begin_tests(src/utils/SiconosTools/test DEPS "CPPUNIT::CPPUNIT")
  new_test_1(SOURCES SiconosGraphTest.cpp ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES SiconosVisitorTest.cpp ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES  SiconosPropertiesTest.cpp ${SIMPLE_TEST_MAIN})

  # ---- Modeling tools ---
  begin_tests(src/modelingTools/test DEPS "numerics;CPPUNIT::CPPUNIT")
  new_test_1(SOURCES FirstOrderNonLinearDSTest.cpp ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES FirstOrderLinearDSTest.cpp ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES FirstOrderLinearTIRTest.cpp ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES FirstOrderLinearRTest.cpp ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES FirstOrderType1RTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES LagrangianLinearTIRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES LagrangianScleronomousRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES LagrangianRheonomousRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES LagrangianCompliantRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES LagrangianCompliantLinearTIRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES LagrangianDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES LagrangianLinearTIDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES NewtonEulerDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test_1(SOURCES NonSmoothDynamicalSystemTest.cpp  ${SIMPLE_TEST_MAIN})
  
  # ---- Simulation tools ---
  begin_tests(src/simulationTools/test DEPS "numerics;CPPUNIT::CPPUNIT")
  new_test_1(SOURCES OSNSPTest.cpp ${SIMPLE_TEST_MAIN})
  if(HAS_FORTRAN)
    new_test_1(SOURCES ZOHTest.cpp ${SIMPLE_TEST_MAIN})
  endif()

 endif()
