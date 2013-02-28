# a ctest driver for nightly builds
# as ctest cannot go through differents projects, this driver is called from Build/CMakeLists.txt
# if -DON_DASHBOARD=1 is passed to cmake

if(NOT MODULE)
  set(MODULE "Numerics")
endif()

if(NOT MODE)
  message(STATUS "set mode to Experimental")
  set(MODE "Experimental")
endif()

set(BUILDNAME ${MODULE})

set(CTEST_SOURCE_DIRECTORY "${SOURCE_DIRECTORY}/${MODULE}")
set(CTEST_BINARY_DIRECTORY "${BINARY_DIRECTORY}")

# here we should have all Siconos sources
set(CMAKE_MODULE_PATH ${SOURCE_DIRECTORY}/${MODULE}/cmake)

if(NOT BUILD_TYPE)
  set(CTEST_BUILD_CONFIGURATION "Profiling")
else()
  set(CTEST_BUILD_CONFIGURATION ${BUILD_TYPE})
endif()

set(CMAKE_BUILD_TYPE ${CTEST_BUILD_CONFIGURATION})
set(WITH_PIPOL_TARGETS)
set(_PROJECT_NAME ${MODULE})
include(Pipol)
include(SiconosCTest)

set(CMAKE_BINARY_DIR ${CTEST_BINARY_DIRECTORY})

if(NOT SITE)
  site_name(SITE)
endif()

message(STATUS "submission site is ${SITE}")
set(CTEST_SITE "${SITE}")

message(STATUS "build name is ${BUILDNAME}")
set(CTEST_BUILD_NAME "${BUILDNAME}")

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
if(FROM_REPO)
  set(CTEST_BUILD_OPTIONS "-DWITH_GIT=ON")
endif()

if(CTEST_BUILD_CONFIGURATION MATCHES "Profiling")
  set(WITH_MEMCHECK TRUE)
  set(WITH_COVERAGE TRUE)
endif()

#######################################################################

#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

find_program(CTEST_GIT_COMMAND NAMES git)
find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)

set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all") 
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/cmake/valgrind.supp)

set(CTEST_NOTES_FILES ${CTEST_BINARY_DIRECTORY}/Testing/Notes/Build)

STRING(REGEX REPLACE "__00__" " " EXTRA_CMAKE_ARGS_L "${EXTRA_CMAKE_ARGS}" )

# source directory should be set by meta CMakeLists.txt
#if(NOT EXISTS "${SOURCE_DIRECTORY}")
#  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone git+ssh://${GFORGE_USER}@scm.gforge.inria.fr//gitroot/siconos/siconos.git ${SOURCE_DIRECTORY}")
#endif()

set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION} -DBUILDNAME:STRING=${MODULE} ${CTEST_BUILD_OPTIONS} ${EXTRA_CMAKE_ARGS_L}")

if(WITH_TESTING)
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTING:BOOL=ON")
  if(WITH_COVERAGE)
    set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTS_COVERAGE:BOOL=ON")
  endif()
endif()
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")

ctest_start("${MODE}")
if(FROM_REPO)
  ctest_update()
endif()
ctest_configure()
ctest_build()
ctest_test()
if (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
  ctest_coverage()
endif (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  ctest_memcheck()
endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
ctest_submit()

