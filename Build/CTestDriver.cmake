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

set(CTEST_SOURCE_DIRECTORY "${SOURCE_DIRECTORY}/${MODULE}")
set(CTEST_BINARY_DIRECTORY "${BINARY_DIRECTORY}")

# here we should have all Siconos sources
set(CMAKE_MODULE_PATH ${SOURCE_DIRECTORY}/${MODULE}/cmake)

set(CTEST_BUILD_CONFIGURATION "Profiling")

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
set(CTEST_BUILD_OPTIONS "-DWITH_GIT=ON")

set(WITH_MEMCHECK TRUE)
set(WITH_COVERAGE TRUE)

#######################################################################

#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

find_program(CTEST_GIT_COMMAND NAMES git)
find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)

set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all") 
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/cmake/valgrind.supp)

set(CTEST_NOTES_FILES "${CMAKE_BINARY_DIR}/Testing/Notes/Build")

# source directory should be set by meta CMakeLists.txt
#if(NOT EXISTS "${SOURCE_DIRECTORY}")
#  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone git+ssh://${GFORGE_USER}@scm.gforge.inria.fr//gitroot/siconos/siconos.git ${SOURCE_DIRECTORY}")
#endif()

set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTING:BOOL=ON -DWITH_TESTS_COVERAGE:BOOL=ON ${CTEST_BUILD_OPTIONS}")
#set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -G\"${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")

ctest_start("${MODE}")
ctest_update()
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

