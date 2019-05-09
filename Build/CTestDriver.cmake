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
set(CTEST_MODULE_PATH ${SOURCE_DIRECTORY}/${MODULE}/cmake)

message(STATUS "MODULE is: ${MODULE}")
message(STATUS "MODE is: ${MODE}")
message(STATUS "SOURCE_DIRECTORY is: ${SOURCE_DIRECTORY}")
message(STATUS "BINARY_DIRECTORY is: ${BINARY_DIRECTORY}")

message(STATUS "cmake module path is: ${CMAKE_MODULE_PATH}")


if(NOT BUILD_TYPE)
  set(CTEST_BUILD_CONFIGURATION "Profiling")
else()
  set(CTEST_BUILD_CONFIGURATION ${BUILD_TYPE})
endif()

message(STATUS "BUILD_TYPE is: ${CTEST_BUILD_CONFIGURATION}")

set(CMAKE_BUILD_TYPE ${CTEST_BUILD_CONFIGURATION})
set(CMAKE_SOURCE_DIR ${SOURCE_DIRECTORY})

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

if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()

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

set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all") 
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/cmake/valgrind.supp)

set(CTEST_NOTES_FILES ${CTEST_BINARY_DIRECTORY}/Testing/Notes/Build)
STRING(REGEX REPLACE "__00__" " " EXTRA_CMAKE_ARGS_L "${EXTRA_CMAKE_ARGS}" )
STRING(REGEX REPLACE "__11__" "\"" EXTRA_CMAKE_ARGS_L "${EXTRA_CMAKE_ARGS_L}" )
STRING(REGEX REPLACE "__22__" ":" EXTRA_CMAKE_ARGS_L "${EXTRA_CMAKE_ARGS_L}" )
STRING(REGEX REPLACE "___" ";" EXTRA_CMAKE_ARGS_L "${EXTRA_CMAKE_ARGS_L}" )

if(TEST_TIMEOUT)
  set(CTEST_TEST_TIMEOUT ${TEST_TIMEOUT})
endif(TEST_TIMEOUT)

#MESSAGE( "EXTRA_CMAKE_ARGS_L :: ${EXTRA_CMAKE_ARGS_L}")


# source directory should be set by meta CMakeLists.txt
#if(NOT EXISTS "${SOURCE_DIRECTORY}")
#  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone git+ssh://${GFORGE_USER}@scm.gforge.inria.fr//gitroot/siconos/siconos.git ${SOURCE_DIRECTORY}")
#endif()

set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\" -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION} -DBUILDNAME:STRING=${MODULE} ${CTEST_BUILD_OPTIONS} ${EXTRA_CMAKE_ARGS_L}")

if(WITH_TESTING)
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTING:BOOL=ON")
  if(WITH_COVERAGE)
    set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTS_COVERAGE:BOOL=ON")
  endif()
endif()
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")


set( dashboard_cache "
     BUILD_TESTING:BOOL=ON
     BUILD_EXAMPLES:BOOL=ON
     BUILD_SHARED_LIBS:BOOL=ON
     MEMORYCHECK_COMMAND:FILEPATH=${CTEST_MEMORYCHECK_COMMAND}
     MEMORYCHECK_COMMAND_OPTIONS:STRING=${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
     MEMORYCHECK_SUPPRESSIONS_FILE:FILEPATH=${CTEST_MEMORYCHECK_SUPPRESSIONS_FILE}
     "
  )

include(ProcessorCount)
ProcessorCount(N)
#if(NOT N EQUAL 0)
  #  IF(NOT MODULE MATCHES "IO")
  #    set(CTEST_BUILD_FLAGS -j${N})
  #  ENDIF()
  #  set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
  #endif()

ctest_start("${MODE}")
if(FROM_REPO)
  ctest_update()
endif()
ctest_configure()
ctest_build()
ctest_test(PARALLEL_LEVEL ${N})
if (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
  ctest_coverage()
endif (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  ctest_memcheck()
endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
ctest_submit()

