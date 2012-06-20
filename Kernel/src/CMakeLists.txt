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
  utils
  utils/boost-bindings
  utils/contract
  utils/SiconosAlgebra
  utils/SiconosException
  utils/SiconosMemory
  utils/SiconosSharedLibrary
  utils/SiconosTools
  xml
  model
  plugin
  modelingTools
  simulationTools
  controlTools
  multiBodyTools
  multiBodyTools/2D
  multiBodyTools/3D
  multiBodyTools/joints
)

#
# Unstable sources
#
# If you want to compile those, don't forget to uncomment the
# include in SimulationTools.hpp
SET(SiconosKernel_Unstable_SRCS 
  simulationTools/PrimalFrictionContact.cpp
  simulationTools/MLCP2.cpp
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
  FILE(GLOB BOOST_HDRS_${_D} utils/boost-bindings/boost/numeric/bindings/${_D}/*.hpp 
    utils/boost-bindings/boost/numeric/bindings/${_D}/*.h)
  INSTALL(FILES ${BOOST_HDRS_${_D}} 
    DESTINATION include/Siconos/Kernel/boost/numeric/bindings/${_D})
ENDFOREACH(_D amos;atlas;blas;cblas;lapack;traits;umfpack)

# config.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/KernelConfig.h DESTINATION include/Siconos/Kernel)

if(CMAKE_CXX_COMPILER)
 include(TestCXXAcceptsFlag)
 CHECK_CXX_ACCEPTS_FLAG("-diag-disable 654" disable_warning)
 CHECK_CXX_ACCEPTS_FLAG("-D__aligned__=ignored" DALIGNEDERROR)
 if(NOT CMAKE_COMPILER_IS_GNUCXX)
  if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
   # Disable warnings with intel compiler due (mainly) to visitors visit function overloading
   if(disable_warning)
    append_cxx_flags("-diag-disable 654")
   endif(disable_warning)
   # Error on intel compiler, see: http://software.intel.com/en-us/forums/showthread.php?t=65041
   # This issue have been solved with ICC >= 12.1
   if(NOT ICCOK)
    if(DALIGNEDERROR)
     append_cxx_flags("-D__aligned__=ignored")
    endif(DALIGNEDERROR)
   endif(NOT ICCOK)
  endif(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
 endif(NOT CMAKE_COMPILER_IS_GNUCXX)
 IF(C_HAVE_WALL)
  APPEND_FLAGS(CMAKE_CXX_FLAGS "-Wall")
 ENDIF(C_HAVE_WALL)
 IF(C_HAVE_WEXTRA)
  APPEND_FLAGS(CMAKE_CXX_FLAGS "-Wextra -Wno-unused-parameter")
 ENDIF(C_HAVE_WEXTRA)
 IF(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
   # stupid warning
  APPEND_FLAGS(CMAKE_CXX_FLAGS "-Wno-string-plus-int")
 ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

endif(CMAKE_CXX_COMPILER)

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
IF(NOT WITHOUT_ATLAS)
COMPILE_WITH(ATLAS)
ENDIF(NOT WITHOUT_ATLAS)
COMPILE_WITH(GMP REQUIRED)

IF(WITH_TESTING)
  COMPILE_WITH(CppUnit REQUIRED)
ENDIF(WITH_TESTING)

# contact detection
IF(WITH_BULLET)
  COMPILE_WITH(Bullet REQUIRED)
  IF(BULLET_FOUND)
    SET(HAVE_BULLET TRUE)
    MESSAGE( STATUS " Bullet include dirs : ${BULLET_INCLUDE_DIRS}" )
    SET(SiconosKernel_DIRS ${SiconosKernel_DIRS} multiBodyTools/contactDetection/bullet)
    IF(BULLET_USE_DOUBLE_PRECISION)
      APPEND_FLAGS(CMAKE_CXX_FLAGS "-DBT_USE_DOUBLE_PRECISION")
    ENDIF(BULLET_USE_DOUBLE_PRECISION)
  ENDIF(BULLET_FOUND)
ENDIF(WITH_BULLET)

# Plugin library for tests
FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/plugin/test)
# MODULE rather than SHARED for Macosx compatibility. 
ADD_LIBRARY(TestPlugin MODULE plugin/test/TestPlugin.cpp)
SET_TARGET_PROPERTIES(TestPlugin 
  PROPERTIES PREFIX ""
  OUTPUT_NAME plugin/test/TestPlugin)
# Note Franck. Is it really necessary to install it? I Don't think so ...
#INSTALL(TARGETS TestPlugin 
#  LIBRARY DESTINATION share/${PROJECT_PACKAGE_NAME})

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
    BlockMatrixTest.cpp  SimpleMatrixTest.cpp BlockVectorTest.cpp  SiconosVectorTest.cpp)
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

  NEW_TEST(testSiconosTools SiconosGraphTest.cpp SiconosVisitorTest.cpp SiconosPropertiesTest.cpp)

  END_TEST()

  BEGIN_TEST(multiBodyTools/test)

  NEW_TEST(testMultiBodyTools MultiBodyTest.cpp)

  END_TEST()

  # Simulation tests
  BEGIN_TEST(simulationTools/test)

  NEW_TEST(testSimulationTools ZOHTest.cpp)

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


# To have a full print of variables, use -DPRINT_ENV=ON
MESSAGE(STATUS "\n============================================ Summary ============================================")
MESSAGE(STATUS "${PROJECT_NAME} version ${VERSION} is now ready for compilation and installation.")
MESSAGE(STATUS "To proceed run 'make' and 'make install' and optionaly 'make test'.")
MESSAGE(STATUS " C++ Compiler : ${CMAKE_CXX_COMPILER}")
MESSAGE(STATUS " C Compiler : ${CMAKE_C_COMPILER}")
MESSAGE(STATUS " Fortran Compiler : ${CMAKE_Fortran_COMPILER}")
MESSAGE(STATUS " Compilation mode is : ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS " Code Sources are in : ${CMAKE_SOURCE_DIR}")
MESSAGE(STATUS " The project will be installed in ${CMAKE_INSTALL_PREFIX}")
MESSAGE(STATUS " The Numerics library used is : ${SiconosNumerics_FOUND}")
MESSAGE(STATUS " To get more information about dependencies, config or else, ")
MESSAGE(STATUS "check CMakeCache.txt file or re-run cmake with -DPRINT_ENV=ON.")
MESSAGE(STATUS "=================================================================================================\n")
