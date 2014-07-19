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
  utils/numeric_bindings/
  utils/contract
  utils/SiconosAlgebra
  utils/SiconosException
  utils/SiconosMemory
  utils/SiconosSharedLibrary
  utils/SiconosTools
  model
  plugin
  modelingTools
  simulationTools
  controlTools
)

#
# Unstable sources
#
# If you want to compile those, don't forget to uncomment the
# include in SimulationTools.hpp
SET(SiconosKernel_Unstable_SRCS 
  simulationTools/GlobalFrictionContact.cpp
  simulationTools/MLCP2.cpp
  )


#
# Installation layout
#
SET(SiconosKernel_INSTALL_INCLUDE_DIR include/Siconos/Kernel)

# Atlas and boost headers
#FILE(GLOB ATLAS_HDRS utils/AtlasLocal/atlas/*.h)
#INSTALL(FILES ${ATLAS_HDRS} 
#  DESTINATION include/Siconos/Kernel/AtlasLocal/atlas)

FOREACH(_D blas;lapack;traits;traits/detail;umfpack;std;ublas;ublas/detail;detail)
  FILE(GLOB BOOST_HDRS_${_D} utils/numeric_bindings/boost/numeric/bindings/${_D}/*.hpp 
    utils/numeric_bindings/boost/numeric/bindings/${_D}/*.h)
  INSTALL(FILES ${BOOST_HDRS_${_D}} 
    DESTINATION include/Siconos/Kernel/boost/numeric/bindings/${_D})
ENDFOREACH()
FILE(GLOB BOOST_HDRS utils/numeric_bindings/boost/numeric/bindings/*.hpp) 

INSTALL(FILES ${BOOST_HDRS} 
  DESTINATION include/Siconos/Kernel/boost/numeric/bindings/)

# config.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/KernelConfig.h DESTINATION include/Siconos/Kernel)

IF(CMAKE_CXX_COMPILER)

 IF(NOT CMAKE_COMPILER_IS_GNUCXX)
  IF(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
   # Disable warnings with intel compiler due (mainly) to visitors visit function overloading
   IF(CXX_HAVE_DIAG_DISABLE_654)
    APPEND_CXX_FLAGS("-diag-disable 654")
   ENDIF(CXX_HAVE_DIAG_DISABLE_654)
   # Error on intel compiler, see: http://software.intel.com/en-us/forums/showthread.php?t=65041
   # This issue have been solved with ICC >= 12.1
   if(NOT ICCOK)
    if(CXX_HAVE_D__ALIGNED__IGNORED)
     APPEND_CXX_FLAGS("-D__aligned__=ignored")
    endif(CXX_HAVE_D__ALIGNED__IGNORED)
   endif(NOT ICCOK)
  endif(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
 endif(NOT CMAKE_COMPILER_IS_GNUCXX)
 # way too verbose with MSVC
 IF(CXX_HAVE_WALL AND NOT MSVC)
  APPEND_CXX_FLAGS("-Wall")
 ENDIF(CXX_HAVE_WALL AND NOT MSVC)
 IF(CXX_HAVE_WEXTRA)
  APPEND_CXX_FLAGS("-Wextra -Wno-unused-parameter")
 ENDIF(CXX_HAVE_WEXTRA)
 IF(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
   # stupid warning
  APPEND_CXX_FLAGS("-Wno-string-plus-int")
 ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

 IF(CXX_HAVE_WOVERLOADED_VIRTUAL)
   APPEND_CXX_FLAGS("-Woverloaded-virtual")
 ENDIF(CXX_HAVE_WOVERLOADED_VIRTUAL)

 IF(CXX_HAVE_IMPL AND NOT MSVC)
  APPEND_CXX_FLAGS("-Werror=implicit-function-declaration")
 ENDIF(CXX_HAVE_IMPL AND NOT MSVC)

 IF(CXX_HAVE_UNREACH AND NOT MSVC)
    APPEND_CXX_FLAGS("-Werror=unreachable-code")
 ENDIF(CXX_HAVE_UNREACH AND NOT MSVC)

 IF(CXX_HAVE_CONV AND NOT MSVC)
   APPEND_CXX_FLAGS("-Werror=conversion -Wno-sign-conversion")
   APPEND_CXX_FLAGS("-Wno-error=sign-conversion -Wno-error=shorten-64-to-32")
 ENDIF(CXX_HAVE_CONV AND NOT MSVC)

 IF(CXX_HAVE_MISS AND NOT MSVC)
    APPEND_CXX_FLAGS("-Wmissing-declarations")
 ENDIF(CXX_HAVE_MISS AND NOT MSVC)

ENDIF(CMAKE_CXX_COMPILER)

INCLUDE(cxxVersion)
detect_cxx_version(CXX_VERSION)


COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(Boost REQUIRED)
# FIND_PACKAGE( ... COMPONENTS <component>) does not mention <component> in failure message!
# a workaround 

#FIND_PACKAGE(Boost COMPONENTS graph)

#IF(NOT Boost_GRAPH_FOUND)
#  MESSAGE(STATUS "WARNING boost graph library not found.")
#ENDIF(NOT Boost_GRAPH_FOUND)

compile_with(BLAS REQUIRED)
compile_with(LAPACK REQUIRED)

COMPILE_WITH(GMP REQUIRED)

IF(WITH_TESTING)
  # We don't use COMPILE_WITH since we don't want to link cppunit with the
  # kernel library
  find_package(CppUnit REQUIRED)
  set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})
  set(TEST_INCLUDE_DIR ${TEST_INCLUDE_DIR} ${CPPUNIT_INCLUDE_DIR})
ENDIF(WITH_TESTING )


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

# for dlopen or equivalent (used for plugins)
SET(SiconosKernel_LINK_LIBRARIES ${SiconosKernel_LINK_LIBRARIES}
 ${CMAKE_DL_LIBS})

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  STRING(REGEX REPLACE "\\\\" "/" ENV_PATH "$ENV{PATH}")
  STRING(REGEX REPLACE "\;" "\\\;" ENV_PATH "${ENV_PATH}")
ENDIF()

LIBRARY_PROJECT_SETUP()

# wrappers for tests are needed
# XXX maybe it's time to get rid of this -- xhub
IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
  SET(TEST_WRAP FALSE)
ELSE()
  SET(TEST_WRAP FALSE)
ENDIF(CROSSCOMPILING_LINUX_TO_WINDOWS)

IF(MSVC)
  SET(CMAKE_CXX_CREATE_SHARED_MODULE
"<CMAKE_LINKER> ${CMAKE_CL_NOLOGO} <OBJECTS> ${CMAKE_START_TEMP_FILE} /out:<TARGET> /pdb:<TARGET_PDB> /dll /version:<TARGET_VERSION_MAJOR>.<TARGET_VERSION_MINOR> <LINK_FLAGS> <LINK_LIBRARIES> ${CMAKE_END_TEMP_FILE}")
  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/siconos/cmake-vc.bat.in
    ${CMAKE_BINARY_DIR}/cmake-vc.bat)
  INSTALL(FILES ${CMAKE_BINARY_DIR}/cmake-vc.bat DESTINATION bin)
  INSTALL(FILES ${CMAKE_SOURCE_DIR}/siconos/siconos.bat DESTINATION bin)
ENDIF(MSVC)
#
# tests with cppunit
#
IF(WITH_TESTING)

  # the main test driver
  SET(TEST_MAIN model/test/TestMain.cpp)

  # For Windows
  SET(PATH_FOR_PLUGIN
   ".\;${CMAKE_CURRENT_BINARY_DIR}/plugin\;${CMAKE_CURRENT_BINARY_DIR}/plugin/test")
  #  GET_FILENAME_COMPONENT(SICONOS_LIB_PATH SiconosNumerics_FOUND PATH)
  SET(SICONOS_LIB_PATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")

  # Siconos Algebra
  BEGIN_TEST(utils/SiconosAlgebra/test)

  NEW_TEST(testSiconosAlgebra
    BlockMatrixTest.cpp  SimpleMatrixTest.cpp BlockVectorTest.cpp  SiconosVectorTest.cpp EigenProblemsTest.cpp AlgebraToolsTest.cpp)
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

  # Global tests
  BEGIN_TEST(global_tests)

  ADD_LIBRARY(gTestPlugin MODULE global_tests/gTestPlugin.cpp)
  SET_TARGET_PROPERTIES(gTestPlugin 
    PROPERTIES PREFIX "" 
    OUTPUT_NAME global_tests/gTestPlugin)
  
  NEW_TEST(testGlobal
    DiodeBridge.cpp FunctionTest.cpp  GlobalTests.cpp)
  
  END_TEST()
  
  BEGIN_TEST(utils/SiconosTools/test)

  NEW_TEST(testSiconosTools SiconosGraphTest.cpp SiconosVisitorTest.cpp SiconosPropertiesTest.cpp)

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
# double expansion required
SET(MY_AT "@")
IF(CROSSCOMPILING_LINUX_TO_WINDOWS)
  SET(EXTRA_EXT ".a")
ELSE()
  SET(EXTRA_EXT)
ENDIF()
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/Settings.cmake.in 
 ${CMAKE_BINARY_DIR}/Settings_tmp.cmake @ONLY)
CONFIGURE_FILE(${CMAKE_BINARY_DIR}/Settings_tmp.cmake
 ${CMAKE_BINARY_DIR}/Settings.cmake @ONLY)

# generated cmake macros
INSTALL(FILES 
  ${CMAKE_BINARY_DIR}/LibraryDependencies.cmake 
  ${CMAKE_BINARY_DIR}/Settings.cmake
  DESTINATION share/${PROJECT_PACKAGE_NAME}/cmake)

# other needed cmake macros
INSTALL(FILES 
  ${CMAKE_SOURCE_DIR}/cmake/FindQGLViewer.cmake
  DESTINATION share/${PROJECT_PACKAGE_NAME}/cmake)

INSTALL(FILES 
  ${CMAKE_SOURCE_DIR}/cmake/FindSiconosNumerics.cmake
  DESTINATION share/${PROJECT_PACKAGE_NAME}/cmake)

INSTALL(FILES 
  ${CMAKE_SOURCE_DIR}/cmake/FindSiconosIO.cmake
  DESTINATION share/${PROJECT_PACKAGE_NAME}/cmake)

INSTALL(FILES 
  ${CMAKE_SOURCE_DIR}/cmake/FindSiconosMechanics.cmake
  DESTINATION share/${PROJECT_PACKAGE_NAME}/cmake)

INSTALL(FILES 
  ${CMAKE_SOURCE_DIR}/cmake/SiconosTools.cmake
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

# print all features
FEATURE_SUMMARY(WHAT ALL)

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
MESSAGE(STATUS " The Numerics library used is : ${SiconosNumerics_LIBRARIES}")
MESSAGE(STATUS " The Numerics headers are located in : ${SiconosNumerics_INCLUDE_DIRS}")
message(STATUS " Blas is ${WITH_BLAS}. Libraries : ${BLAS_LIBRARIES}. Header : ${BLAS_INCLUDE_DIRS}/${BLAS_HEADER}.")
message(STATUS " Lapack is ${WITH_LAPACK}. Libraries : ${LAPACK_LIBRARIES}. Header : ${LAPACK_INCLUDE_DIRS}/${LAPACK_HEADER}.")
MESSAGE(STATUS " To get more information about dependencies, config or else, ")
MESSAGE(STATUS "check CMakeCache.txt file or re-run cmake with -DPRINT_ENV=ON.")
MESSAGE(STATUS "=================================================================================================\n")

WRITE_NOTES()


