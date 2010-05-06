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
  plugin
  modelingTools
  simulationTools
  controlTools
  multiBodyTools
  multiBodyTools/2D
  multiBodyTools/3D
  multiBodyTools/joints
  utils
  utils/boost-bindings
  utils/contract
  utils/SiconosAlgebra
  utils/SiconosException
  utils/SiconosMemory
  utils/SiconosSharedLibrary
  utils/SiconosTools
  xml
)

#
# Unstable sources
#
SET(SiconosKernel_Unstable_SRCS 
  simulationTools/PrimalFrictionContact.cpp
  simulationTools/MLCP2.cpp)


#
# Installation layout
#
SET(SiconosKernel_INSTALL_INCLUDE_DIR include/Siconos/Kernel)

# Atlas and boost headers
FILE(GLOB ATLAS_HDRS utils/AtlasLocal/atlas/*.h)
INSTALL(FILES ${ATLAS_HDRS} 
  DESTINATION include/Siconos/Kernel/AtlasLocal/atlas)

FOREACH(_D amos;atlas;blas;cblas;lapack;traits;traits/detail;umfpack)
  FILE(GLOB BOOST_HDRS_${_D} utils/boost-bindings/boost/numeric/bindings/${_D}/*.hpp 
    utils/boost-bindings/boost/numeric/bindings/${_D}/*.h)
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

COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(Boost REQUIRED)
# FIND_PACKAGE( ... COMPONENTS <component>) does not mention <component> in failure message!
# a workaround 

FIND_PACKAGE(Boost COMPONENTS graph)

IF(NOT Boost_GRAPH_FOUND)
  MESSAGE(STATUS "WARNING boost graph library not found.")
ENDIF(NOT Boost_GRAPH_FOUND)

COMPILE_WITH(LibXml2 REQUIRED)
COMPILE_WITH(LAPACK REQUIRED)
COMPILE_WITH(BLAS REQUIRED)
COMPILE_WITH(ATLAS)
COMPILE_WITH(GMP REQUIRED)
COMPILE_WITH(CppUnit)

IF(NOT CPPUNIT_FOUND)
  SET(WITH_TESTING FALSE)
ENDIF(NOT CPPUNIT_FOUND)

FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/plugin/test)

ADD_LIBRARY(TestPlugin MODULE plugin/test/TestPlugin.cpp)
SET_TARGET_PROPERTIES(TestPlugin 
  PROPERTIES PREFIX ""
  OUTPUT_NAME plugin/test/TestPlugin)

INSTALL(TARGETS TestPlugin 
  LIBRARY DESTINATION share/${PROJECT_PACKAGE_NAME})

SET(SiconosKernel_LINKER_LANGUAGE CXX)


# do not skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

LIBRARY_PROJECT_SETUP()

# wrappers for tests are needed
SET(TEST_WRAP TRUE)

#
# tests with cppunit
#
IF(WITH_TESTING)

  # the main test driver
  SET(TEST_MAIN model/test/TestMain.cpp)

  # used in ldwrap.c.in
  SET(LD_LIBRARY_PATH ".:${CMAKE_CURRENT_BINARY_DIR}/plugin:${CMAKE_CURRENT_BINARY_DIR}/plugin/test")

  # Siconos Algebra
  BEGIN_TEST(utils/SiconosAlgebra/test)

  NEW_TEST(testSiconosAlgebra
    BlockMatrixTest.cpp  SimpleMatrixTest.cpp BlockVectorTest.cpp  SimpleVectorTest.cpp)
  END_TEST()
  
  # Siconos Memory
  BEGIN_TEST(utils/SiconosMemory/test)
  
  NEW_TEST(testSiconosMemory SiconosMemoryTest.cpp)
  END_TEST()
  
  # model
  BEGIN_TEST(model/test)
  
  NEW_TEST(testModel ModelTest.cpp)
  
  END_TEST()
  # modeling tools 
  BEGIN_TEST(modelingTools/test)
  
  NEW_TEST(testModelingTools FirstOrderLinearDSTest.cpp FirstOrderLinearTIRTest.cpp FirstOrderLinearRTest.cpp 
    FirstOrderType1RTest.cpp 
    LagrangianLinearTIRTest.cpp LagrangianScleronomousRTest.cpp LagrangianRheonomousRTest.cpp 
    LagrangianCompliantRTest.cpp)
  END_TEST()
  #FirstOrderNonLinearDSTest.cpp FirstOrderLinearDSTest.cpp 
  #LagrangianDSTest.cpp LagrangianLinearTIDSTest.cpp TestMain.cpp)
    
    
  
  # xml
  BEGIN_TEST(xml/test)
  NEW_TEST(testXMLTools )
  
  END_TEST()
  
  # Global tests
  BEGIN_TEST(global_tests)
  
  ADD_LIBRARY(gTestPlugin MODULE global_tests/gTestPlugin.cpp)
  SET_TARGET_PROPERTIES(gTestPlugin 
    PROPERTIES PREFIX "" 
    OUTPUT_NAME global_tests/gTestPlugin)
  
  NEW_TEST(testGlobal
    DiodeBridge.cpp BallBowl.cpp BouncingBall.cpp FunctionTest.cpp  GlobalTests.cpp)
  
  END_TEST()
  
  BEGIN_TEST(utils/SiconosTools/test)

  NEW_TEST(testSiconosGraph SiconosGraphTest.cpp)
  NEW_TEST(testSiconosVisitor SiconosVisitorTest.cpp)


  END_TEST()

ENDIF(WITH_TESTING)

#
# Export
#

# standard export macros
EXPORT_LIBRARY_DEPENDENCIES(LibraryDependencies.cmake)

# installation of exported files
CHECK_INSTALL_INCLUDE_DIRECTORIES()

# specific settings
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/Settings.cmake.in 
  ${CMAKE_BINARY_DIR}/Settings.cmake)

# generated cmake macros
INSTALL(FILES 
  ${CMAKE_BINARY_DIR}/LibraryDependencies.cmake 
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

IF(PRINT_VAR)
  MESSAGE(STATUS)
  PRINT_VAR(${PRINT_VAR})
  MESSAGE(STATUS)
ENDIF(PRINT_VAR)
