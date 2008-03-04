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
  plugin
  simulationTools
  utils
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

LIBRARY_PROJECT_SETUP()

#
# tests
#

# model
TEST_DIRECTORY(model/test)
SET_TEST(test1 ModelTest.cpp TestMain.cpp)

# modeling tools 
TEST_DIRECTORY(modelingTools/test)

SET_TEST(test1 
  FirstOrderType1RTest.cpp FirstOrderLinearRTest.cpp FirstOrderLinearTIRTest.cpp 
  LagrangianScleronomousRTest.cpp LagrangianRheonomousRTest.cpp LagrangianCompliantRTest.cpp 
  LagrangianLinearRTest.cpp FirstOrderNonLinearDSTest.cpp FirstOrderLinearDSTest.cpp 
  LagrangianDSTest.cpp LagrangianLinearTIDSTest.cpp TestMain.cpp)

# xml
TEST_DIRECTORY(xml/test)
SET_TEST(test1 
  TestMain.cpp)

IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)

