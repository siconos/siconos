include(tools4tests)

set(TEST_WRAP FALSE)

if(WITH_${COMPONENT}_TESTING)
  
  # We don't use COMPILE_WITH since we don't want to link cppunit with the
  # kernel library
  find_package(CppUnit REQUIRED)
  set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})
  set(TEST_INCLUDE_DIR ${TEST_INCLUDE_DIR} ${CPPUNIT_INCLUDE_DIR})

  # Plugin library for tests
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/plugin/test)
  # MODULE rather than SHARED for Macosx compatibility. 
  add_library(TestPlugin MODULE src/plugin/test/TestPlugin.cpp)
  set_target_properties(TestPlugin 
    PROPERTIES PREFIX ""
    OUTPUT_NAME TestPlugin)

  # the main test driver
  SET(TEST_MAIN src/test/TestMain.cpp)

  # For Windows
  SET(PATH_FOR_PLUGIN
   ".\;${CMAKE_CURRENT_BINARY_DIR}/src/plugin\;${CMAKE_CURRENT_BINARY_DIR}/src/plugin/test")

 # Siconos Algebra
  BEGIN_TEST(src/utils/SiconosAlgebra/test)

  NEW_TEST(testSiconosAlgebra
    BlockMatrixTest.cpp  SimpleMatrixTest.cpp BlockVectorTest.cpp  SiconosVectorTest.cpp EigenProblemsTest.cpp AlgebraToolsTest.cpp)
  END_TEST()
  
  # Siconos Memory
  BEGIN_TEST(src/utils/SiconosMemory/test)
  
  NEW_TEST(testSiconosMemory SiconosMemoryTest.cpp)
  END_TEST()
  
  # modeling tools 
  BEGIN_TEST(src/modelingTools/test)
  
  NEW_TEST(testModelingTools
    FirstOrderNonLinearDSTest.cpp
    FirstOrderLinearDSTest.cpp
    FirstOrderLinearTIRTest.cpp
    FirstOrderLinearRTest.cpp
    FirstOrderType1RTest.cpp 
    LagrangianLinearTIRTest.cpp
    LagrangianScleronomousRTest.cpp
    LagrangianRheonomousRTest.cpp
    LagrangianCompliantRTest.cpp
    LagrangianCompliantLinearTIRTest.cpp
    LagrangianDSTest.cpp
    LagrangianLinearTIDSTest.cpp
    NewtonEulerDSTest.cpp
    NonSmoothDynamicalSystemTest.cpp)
  END_TEST()
  #FirstOrderNonLinearDSTest.cpp FirstOrderLinearDSTest.cpp 
  #LagrangianDSTest.cpp LagrangianLinearTIDSTest.cpp TestMain.cpp)

  
  BEGIN_TEST(src/utils/SiconosTools/test)

  NEW_TEST(testSiconosTools SiconosGraphTest.cpp SiconosVisitorTest.cpp SiconosPropertiesTest.cpp)

  END_TEST()

  # Simulation tests
  BEGIN_TEST(src/simulationTools/test)

  IF(HAS_FORTRAN)
   NEW_TEST(testSimulationTools OSNSPTest.cpp EulerMoreauTest.cpp LsodarTest.cpp ZOHTest.cpp)
   ELSE()
    NEW_TEST(testSimulationTools OSNSPTest.cpp EulerMoreauTest.cpp)
  ENDIF()
  
  END_TEST()

 endif()
