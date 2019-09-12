include(tools4tests)

if(WITH_${COMPONENT}_TESTING)

  # ---- Siconos Algebra tests ----
  begin_tests(src/utils/SiconosAlgebra/test)

  new_test(
    NAME testSiconosAlgebra
    SOURCES BlockMatrixTest.cpp  SimpleMatrixTest.cpp BlockVectorTest.cpp  SiconosVectorTest.cpp EigenProblemsTest.cpp AlgebraToolsTest.cpp ${SIMPLE_TEST_MAIN}
    DEPS "numerics;CPPUNIT::CPPUNIT;externals"
    )

  # ---- Siconos Memory tests ----
  begin_tests(src/utils/SiconosMemory/test)

  new_test(
    NAME testSiconosMemory
    SOURCES SiconosMemoryTest.cpp ${SIMPLE_TEST_MAIN}
    DEPS "numerics;CPPUNIT::CPPUNIT"
    )

  # ---- Siconos tools tests ----
  begin_tests(src/utils/SiconosTools/test DEPS "CPPUNIT::CPPUNIT")
  new_test(SOURCES SiconosGraphTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES SiconosVisitorTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES  SiconosPropertiesTest.cpp ${SIMPLE_TEST_MAIN})

  # ---- Modeling tools ---
  begin_tests(src/modelingTools/test DEPS "numerics;CPPUNIT::CPPUNIT")
  new_test(SOURCES FirstOrderNonLinearDSTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES FirstOrderLinearDSTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES FirstOrderLinearTIRTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES FirstOrderLinearRTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES FirstOrderType1RTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianLinearTIRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianScleronomousRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianRheonomousRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianCompliantRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianCompliantLinearTIRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianLinearTIDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES NewtonEulerDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES NonSmoothDynamicalSystemTest.cpp  ${SIMPLE_TEST_MAIN})
  
  # ---- Simulation tools ---
  begin_tests(src/simulationTools/test DEPS "numerics;CPPUNIT::CPPUNIT")
  new_test(SOURCES OSNSPTest.cpp ${SIMPLE_TEST_MAIN} DEPS LAPACK::LAPACK)
  if(HAS_FORTRAN)
    new_test(SOURCES ZOHTest.cpp ${SIMPLE_TEST_MAIN} DEPS LAPACK::LAPACK)
  endif()

 endif()
