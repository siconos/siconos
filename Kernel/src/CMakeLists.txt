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
  model
  modelingTools
  simulationTools
  controlTools
  utils
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
# Drawstuff sources
#
SET(SiconosKernel_Drawstuff_SRCS 
  multiBodyTools
  multiBodyTools/drawstuff
  multiBodyTools/textures
  )

#
# Installation layout
#
SET(SiconosKernel_INSTALL_INCLUDE_DIR include/Siconos/Kernel)

# Atlas and boost headers
FILE(GLOB ATLAS_HDRS utils/AtlasLocal/atlas/*.h)
INSTALL(FILES ${ATLAS_HDRS} 
  DESTINATION include/Siconos/Kernel/AtlasLocal/atlas)
FOREACH(_D amos;atlas;blas;cblas;lapack;traits;traits/detail;umfpack)
  FILE(GLOB BOOST_HDRS_${_D} utils/boost/numeric/bindings/${_D}/*.hpp)
  INSTALL(FILES ${BOOST_HDRS_${_D}} 
    DESTINATION include/Siconos/Kernel/boost/numeric/bindings/${_D})
ENDFOREACH(_D amos;atlas;blas;cblas;lapack;traits;umfpack)

# config.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/KernelConfig.h DESTINATION include/Siconos/Kernel)

#
# compilation setup 
#
IF(CXX_HAVE_FRIEND_INJECTION)
  APPEND_FLAGS(CMAKE_CXX_FLAGS "-ffriend-injection")
ENDIF(CXX_HAVE_FRIEND_INJECTION)

# Boost NDEBUG
ADD_DEFINITIONS(-DNDEBUG)
COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(Boost REQUIRED)
COMPILE_WITH(LibXml2 REQUIRED)
COMPILE_WITH(LAPACK REQUIRED)
COMPILE_WITH(BLAS REQUIRED)
COMPILE_WITH(ATLAS)
COMPILE_WITH(GMP REQUIRED)
COMPILE_WITH(CppUnit)
COMPILE_WITH(OpenGL)

IF(OPENGL_FOUND)
  SET(SiconosKernel_DIRS ${SiconosKernel_DIRS} ${SiconosKernel_Drawstuff_SRCS})
ENDIF(OPENGL_FOUND)

FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/plugin/test)
ADD_LIBRARY(DefaultPlugin MODULE plugin/DefaultPlugin.cpp)
SET_TARGET_PROPERTIES(DefaultPlugin 
  PROPERTIES PREFIX ""
  OUTPUT_NAME plugin/DefaultPlugin)

ADD_LIBRARY(TestPlugin MODULE plugin/test/TestPlugin.cpp)
SET_TARGET_PROPERTIES(TestPlugin 
  PROPERTIES PREFIX ""
  OUTPUT_NAME plugin/test/TestPlugin)

INSTALL(TARGETS DefaultPlugin TestPlugin 
  LIBRARY DESTINATION share/${PROJECT_PACKAGE_NAME})

SET(SiconosKernel_LINKER_LANGUAGE CXX)

LIBRARY_PROJECT_SETUP()

#
# tests with cppunit
#
IF(CPPUNIT_FOUND)
  # used in ldwrap.c.in
  SET(LD_LIBRARY_PATH ".:${CMAKE_CURRENT_BINARY_DIR}/plugin:${CMAKE_CURRENT_BINARY_DIR}/plugin/test")

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
  
  ADD_LIBRARY(gTestPlugin MODULE global_tests/gTestPlugin.cpp)
  SET_TARGET_PROPERTIES(gTestPlugin 
    PROPERTIES PREFIX "" 
    OUTPUT_NAME global_tests/gTestPlugin)
  
  NEW_TEST(testGlobal
    BallBowl.cpp  BouncingBall.cpp  DiodeBridge.cpp  FunctionTest.cpp  GlobalTests.cpp  MainTest.cpp)
  
  END_TEST()
  
ENDIF(CPPUNIT_FOUND)

#
# Export
#

# standard export macros
EXPORT_LIBRARY_DEPENDENCIES(LibraryDependencies.cmake)
CMAKE_EXPORT_BUILD_SETTINGS(BuildSettings.cmake)

# installation of exported files
CHECK_INSTALL_INCLUDE_DIRECTORIES()

# specific settings
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/Settings.cmake.in 
  ${CMAKE_BINARY_DIR}/Settings.cmake)

# generated cmake macros
INSTALL(FILES 
  ${CMAKE_BINARY_DIR}/LibraryDependencies.cmake 
  ${CMAKE_BINARY_DIR}/BuildSettings.cmake 
  ${CMAKE_BINARY_DIR}/Settings.cmake
  DESTINATION share/${PROJECT_PACKAGE_NAME}/cmake)

# other needed cmake macros
INSTALL(FILES 
  ${CMAKE_SOURCE_DIR}/cmake/FindQGLViewer.cmake
  DESTINATION share/${PROJECT_PACKAGE_NAME}/cmake)

#
# Log 
#

IF(PRINT_ENV)
  INCLUDE(CMakeLog)
ENDIF(PRINT_ENV)
