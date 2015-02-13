#
# Here comes Control
#
PROJECT(SiconosControl CXX)

#
# VERSION defined in top-level CMakeLists.txt
#
SET(SiconosControl_VERSION "${VERSION}")

#
# Published headers
#
# if <PROJECT>_HDRS unset -> all headers under <PROJECT>_DIRS

#
# Project directories 
#
SET(SiconosControl_DIRS
  .
  utils
  Controller
  Observer
  Sensor
  Simulation
)

#
# Unstable sources
#
SET(SiconosControl_Unstable_SRCS)


#
# Installation layout
#
SET(SiconosControl_INSTALL_INCLUDE_DIR include/Siconos/Control)

# config.h
INSTALL(FILES ${CMAKE_BINARY_DIR}/ControlConfig.h DESTINATION include/Siconos/Control)

INCLUDE(SiconosCompilerFlags)

INCLUDE(cxxVersion)
detect_cxx_version(CXX_VERSION)


COMPILE_WITH(SiconosNumerics REQUIRED)
COMPILE_WITH(SiconosKernel REQUIRED)
COMPILE_WITH(Boost REQUIRED)
# FIND_PACKAGE( ... COMPONENTS <component>) does not mention <component> in failure message!
# a workaround 

IF(WITH_TESTING)
  # We don't use COMPILE_WITH since we don't want to link cppunit with the
  # kernel library
  find_package(CppUnit REQUIRED)
  set(TEST_LIBS ${TEST_LIBS} ${CPPUNIT_LIBRARIES})
  set(TEST_INCLUDE_DIR ${TEST_INCLUDE_DIR} ${CPPUNIT_INCLUDE_DIR})
ENDIF(WITH_TESTING )


# Plugin library for tests
FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/plugin/test)

SET(SiconosControl_LINKER_LANGUAGE CXX)

# for dlopen or equivalent (used for plugins)
SET(SiconosControl_LINK_LIBRARIES ${SiconosControl_LINK_LIBRARIES}
 ${CMAKE_DL_LIBS})

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  STRING(REGEX REPLACE "\\\\" "/" ENV_PATH "$ENV{PATH}")
  STRING(REGEX REPLACE "\;" "\\\;" ENV_PATH "${ENV_PATH}")
ENDIF()

IF(MSVC)
  SET(CMAKE_CXX_CREATE_SHARED_MODULE
"<CMAKE_LINKER> ${CMAKE_CL_NOLOGO} <OBJECTS> ${CMAKE_START_TEMP_FILE} /out:<TARGET> /pdb:<TARGET_PDB> /dll /version:<TARGET_VERSION_MAJOR>.<TARGET_VERSION_MINOR> <LINK_FLAGS> <LINK_LIBRARIES> ${CMAKE_END_TEMP_FILE}")
  # we need CMAKE_NM for the EXPORTS hack
  # On other Module this variable is set through FindBinutils since we ask for
  # a Fortran compiler --xhub
  find_program(CMAKE_NM NAMES ${_CMAKE_TOOLCHAIN_PREFIX}nm HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
ENDIF(MSVC)

LIBRARY_PROJECT_SETUP()

SET(TEST_WRAP FALSE)

#
# tests with cppunit
#
IF(WITH_TESTING)

  # the main test driver
  SET(TEST_MAIN tests/TestMain.cpp)

  # For Windows
  SET(PATH_FOR_PLUGIN
   ".\;${CMAKE_CURRENT_BINARY_DIR}/plugin\;${CMAKE_CURRENT_BINARY_DIR}/plugin/test")
  #  GET_FILENAME_COMPONENT(SICONOS_LIB_PATH SiconosNumerics_FOUND PATH)
  SET(SICONOS_LIB_PATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")

  # Simulation tests
  BEGIN_TEST(tests)

  NEW_TEST(tests PIDTest.cpp SMCTest.cpp)

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
MESSAGE(STATUS " The Control library used is : ${SiconosKernel_LIBRARIES}")
MESSAGE(STATUS " The Control headers are located in : ${SiconosKernel_INCLUDE_DIRS}")
MESSAGE(STATUS " To get more information about dependencies, config or else, ")
MESSAGE(STATUS "check CMakeCache.txt file or re-run cmake with -DPRINT_ENV=ON.")
MESSAGE(STATUS "=================================================================================================\n")

WRITE_NOTES()


