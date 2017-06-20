# a ctest driver for Experimental,Continuous,Nightly builds

set(CMAKE_MODULE_PATH ${CTEST_SOURCE_DIRECTORY}/CI/cmake;${CTEST_SOURCE_DIRECTORY}/CI/config;${CTEST_SOURCE_DIRECTORY}/CI)

include(Tools)

if(NOT MODE)
  message(STATUS "set mode to Experimental")
  set(MODE "Experimental")
endif()

# -- Get config --
# i.e. set extra options/values (cmake -Doption=value ...)
# either from file default.cmake or
# from file CI_CONFIG.cmake 
# --> may set SICONOS_CMAKE_OPTIONS
# --> may set DSICONOS_COMPONENTS
if(CI_CONFIG)
  string(REPLACE "," ";" CI_CONFIG_LIST ${CI_CONFIG})
  foreach(_CI ${CI_CONFIG_LIST})
    include(${_CI})
  endforeach(_CI CI_CONFIG_LIST)
else()
  set(CI_CONFIG default)
  include(${CI_CONFIG})
endif()

string(REPLACE "," "-" CI_CONFIG_NAME ${CI_CONFIG})

foreach(option ${SICONOS_CMAKE_OPTIONS})
  set(CI_CONFIGURE_OPTIONS "${CI_CONFIGURE_OPTIONS} ${option}")
endforeach()

if(NOT CMAKE_WRAPPER)
  set(CMAKE_WRAPPER "cmake")
endif()

set(BUILDNAME Siconos-${CI_CONFIG_NAME})

if(NOT CTEST_SOURCE_DIRECTORY)
  # assume build directory is under source directory
  set(CTEST_SOURCE_DIRECTORY ..)
endif()

if(NOT CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY .)
endif()

set(CMAKE_MODULE_PATH ${CTEST_SOURCE_DIRECTORY}/cmake)

if(NOT BUILD_CONFIGURATION)
  set(CTEST_BUILD_CONFIGURATION "Profiling")
else()
  set(CTEST_BUILD_CONFIGURATION ${BUILD_CONFIGURATION})
endif()

message(STATUS "Siconos CTest driver")
message(STATUS "MODE is: ${MODE}")
message(STATUS "CTEST_SOURCE_DIRECTORY is: ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY is: ${CTEST_BINARY_DIRECTORY}")
message(STATUS "CTEST_MODULE_PATH is: ${CTEST_MODULE_PATH}")
message(STATUS "CTEST_BUILD_CONFIGURATION is: ${CTEST_BUILD_CONFIGURATION}")

set(CMAKE_BUILD_TYPE ${CTEST_BUILD_CONFIGURATION})

set(CMAKE_SOURCE_DIR ${CTEST_SOURCE_DIRECTORY})

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
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY}/)

# !!
#file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY})

find_program(CTEST_GIT_COMMAND NAMES git)
find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)

set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-leak-kinds=definite,possible --track-origins=yes --error-limit=no --gen-suppressions=all") 
#set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all") 
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/cmake/valgrind.supp)

set(CTEST_NOTES_FILES ${CTEST_BINARY_DIRECTORY}/Testing/Notes/Build)

if(TEST_TIMEOUT)
  set(CTEST_TEST_TIMEOUT ${TEST_TIMEOUT})
endif(TEST_TIMEOUT)

set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

set(CTEST_CONFIGURE_COMMAND "${CMAKE_WRAPPER} \"-G${CTEST_CMAKE_GENERATOR}\" -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION} ${CTEST_BUILD_OPTIONS}")

if(WITH_COVERAGE)
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTS_COVERAGE:BOOL=ON")
endif()


set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CI_CONFIGURE_OPTIONS}")
if(SICONOS_COMPONENTS)
  message("ooi ${SICONOS_COMPONENTS}")
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DCOMPONENTS=${SICONOS_COMPONENTS}")
endif()

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
ProcessorCount(NP)
set(CTEST_PARALLEL_LEVEL ${NP})

ctest_start("${MODE}")
if(FROM_REPO)
  ctest_update()
endif()
ctest_configure()
ctest_build()
ctest_test(PARALLEL_LEVEL ${NP})
if (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
  ctest_coverage()
endif (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  ctest_memcheck()
endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
ctest_submit()

