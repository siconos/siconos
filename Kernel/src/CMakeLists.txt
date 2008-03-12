#
# Here comes Kernel
#
PROJECT(SiconosKernel C CXX)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosKernel_VERSION "${VERSION}")

#
# Published headers
#
# if <PROJECT>_HDRS unset -> all headers under <PROJECT>_DIRS

#
# Project directories 
#
SET(SiconosKernel_DIRS
  .
  controlTools
  model
  modelingTools
  simulationTools
  utils
  utils/AtlasLocal/atlas
  utils/boost/numeric/bindings/atlas
  utils/contract
  utils/SiconosAlgebra
  utils/SiconosException
  utils/SiconosMemory
  utils/SiconosSharedLibrary
  utils/SiconosTools
  xml)

#
# Unstable sources
#
SET(SiconosKernel_Unstable_SRCS )

#
# Installation layout
#
SET(SiconosKernel_INSTALL_INCLUDE_DIR include/Siconos/Kernel)

#
# compilation setup 
#
IF(CXX_HAVE_FRIEND_INJECTION)
  APPEND_FLAGS(CMAKE_CXX_FLAGS "-ffriend-injection")
ENDIF(CXX_HAVE_FRIEND_INJECTION)

COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(Boost REQUIRED)
COMPILE_WITH(LibXml2 REQUIRED)
COMPILE_WITH(BLAS REQUIRED)
COMPILE_WITH(LAPACK REQUIRED)
COMPILE_WITH(ATLAS REQUIRED)
COMPILE_WITH(GMP REQUIRED)
COMPILE_WITH(CppUnit)

FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/plugin/test)
ADD_LIBRARY(DefaultPlugin SHARED plugin/DefaultPlugin.cpp)
SET_TARGET_PROPERTIES(DefaultPlugin 
  PROPERTIES PREFIX ""
  LINK_FLAGS -module
  OUTPUT_NAME plugin/DefaultPlugin)

ADD_LIBRARY(TestPlugin SHARED plugin/test/TestPlugin.cpp)
SET_TARGET_PROPERTIES(TestPlugin 
  PROPERTIES PREFIX ""
  LINK_FLAGS -module
  OUTPUT_NAME plugin/test/TestPlugin)

INSTALL(TARGETS DefaultPlugin TestPlugin DESTINATION share/${PROJECT_PACKAGE_NAME})

LIBRARY_PROJECT_SETUP()


#
# tests
#

# Siconos Algebra
BEGIN_TEST(utils/SiconosAlgebra/test)

NEW_TEST(testSiconosAlgebra
  BlockMatrixTest.cpp  SimpleMatrixTest.cpp BlockVectorTest.cpp  SimpleVectorTest.cpp TestMain.cpp)

END_TEST()

# Siconos Memory
BEGIN_TEST(utils/SiconosMemory/test)

NEW_TEST(testSiconosMemory
  SiconosMemoryTest.cpp TestMain.cpp)

END_TEST()

# model
BEGIN_TEST(model/test)

NEW_TEST(testModel ModelTest.cpp TestMain.cpp)

END_TEST()
# modeling tools 
BEGIN_TEST(modelingTools/test)

NEW_TEST(testModelingTools
  FirstOrderType1RTest.cpp FirstOrderLinearRTest.cpp FirstOrderLinearTIRTest.cpp 
  LagrangianScleronomousRTest.cpp LagrangianRheonomousRTest.cpp LagrangianCompliantRTest.cpp 
  LagrangianLinearRTest.cpp FirstOrderNonLinearDSTest.cpp FirstOrderLinearDSTest.cpp 
  LagrangianDSTest.cpp LagrangianLinearTIDSTest.cpp TestMain.cpp)

END_TEST()


# xml
BEGIN_TEST(xml/test)
NEW_TEST(testXMLTools
    TestMain.cpp)

END_TEST()

# Global tests
BEGIN_TEST(global_tests)

ADD_LIBRARY(gTestPlugin SHARED global_tests/gTestPlugin.cpp)
SET_TARGET_PROPERTIES(gTestPlugin 
  PROPERTIES PREFIX "" 
  LINK_FLAGS -module
  OUTPUT_NAME global_tests/gTestPlugin)

NEW_TEST(testGlobal
  BallBowl.cpp  BouncingBall.cpp  DiodeBridge.cpp  FunctionTest.cpp  GlobalTests.cpp  MainTest.cpp)

END_TEST()


#
#
#

IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)
