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

NEW_TEST(test1 ModelTest.cpp TestMain.cpp)

END_TEST()
# modeling tools 
BEGIN_TEST(modelingTools/test)

NEW_TEST(test1 
  FirstOrderType1RTest.cpp FirstOrderLinearRTest.cpp FirstOrderLinearTIRTest.cpp 
  LagrangianScleronomousRTest.cpp LagrangianRheonomousRTest.cpp LagrangianCompliantRTest.cpp 
  LagrangianLinearRTest.cpp FirstOrderNonLinearDSTest.cpp FirstOrderLinearDSTest.cpp 
  LagrangianDSTest.cpp LagrangianLinearTIDSTest.cpp TestMain.cpp)

END_TEST()


# xml
BEGIN_TEST(xml/test)
NEW_TEST(test1 
    TestMain.cpp)

END_TEST()

IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)
