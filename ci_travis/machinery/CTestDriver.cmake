# a ctest driver for Experimental,Continuous,Nightly builds

set(CMAKE_MODULE_PATH ${CTEST_SOURCE_DIRECTORY}/CI/cmake;${CTEST_SOURCE_DIRECTORY}/CI/config;${CTEST_SOURCE_DIRECTORY}/CI)

include(Tools)

# -- Get config --
# i.e. set extra options/values (cmake -Doption=value ...)
# either from file default.cmake or
# from file CI_CONFIG.cmake
# --> may set SICONOS_CMAKE_OPTIONS
# --> may set DSICONOS_COMPONENTS
if(CI_CONFIG)
  string(REPLACE "," ";" CI_CONFIG_LIST ${CI_CONFIG})
  foreach(_CI IN LISTS CI_CONFIG_LIST)
    include(${_CI})
  endforeach()
else()
  set(CI_CONFIG default)
  include(${CI_CONFIG})
endif()

string(REPLACE "," "-" CI_CONFIG_NAME ${CI_CONFIG})

set(BUILDNAME Siconos-${CI_CONFIG_NAME})

if(NOT CTEST_SOURCE_DIRECTORY)
  # assume build directory is under source directory
  set(CTEST_SOURCE_DIRECTORY ..)
endif()

# --- Configure setup ---
# - Top level build directory -
# If not specified : current dir.
if(NOT CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY .)
endif()

# Current testing model. Priority: 
# Nightly
# Continuous -> set in task (ci_task.py)
# Experimental : default
if(NOT MODE)
  message(STATUS "set mode to Experimental")
  set(MODE "Experimental")
endif()

if(NOT CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY .)
endif()

if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()

# set(CMAKE_MODULE_PATH ${CTEST_SOURCE_DIRECTORY}/cmake)

include(ProcessorCount)
ProcessorCount(NP)
set(CTEST_PARALLEL_LEVEL ${NP})
if(NOT NP EQUAL 0)
set(CTEST_BUILD_FLAGS -j${NP})
set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${NP})
endif()

set(CMAKE_BUILD_TYPE ${CTEST_BUILD_CONFIGURATION})

set(CMAKE_SOURCE_DIR ${CTEST_SOURCE_DIRECTORY})

# include(SiconosCTest)

set(CMAKE_BINARY_DIR ${CTEST_BINARY_DIRECTORY})

if(NOT SITE)
  site_name(SITE)
endif()

message(STATUS "submission site is ${SITE}")
set(CTEST_SITE "${SITE}")

message(STATUS "build name is ${BUILDNAME}")
set(CTEST_BUILD_NAME "${BUILDNAME}")


if(FROM_REPO)
  set(CTEST_BUILD_OPTIONS "-DWITH_GIT=ON")
endif()

if(CTEST_BUILD_CONFIGURATION MATCHES "Profiling")
  set(WITH_MEMCHECK TRUE)
  set(WITH_COVERAGE TRUE)
endif()

#######################################################################
# this usually fails for some reasons and ctest may returns a fail code.
# ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY}/)
# cf discussions here:
# https://gitlab.kitware.com/cmake/cmake/issues/17000

if(SUBMIT EQUAL 0)
  if(CTEST_BINARY_DIRECTORY)
    file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY})
  endif()
endif()

find_program(CTEST_GIT_COMMAND NAMES git)
find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)

set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-leak-kinds=definite,possible --track-origins=yes --error-limit=no --gen-suppressions=all")
#set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all")
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/cmake/valgrind.supp)

#set(CTEST_NOTES_FILES ${CTEST_BINARY_DIRECTORY}/Testing/Notes/Build)

if(TEST_TIMEOUT)
  set(CTEST_TEST_TIMEOUT ${TEST_TIMEOUT})
endif()

set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")


ctest_start("${MODE}")

# Set CTEST_CONFIGURE_COMMAND to cmake followed by siconos options 
set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\" -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION} ${CTEST_BUILD_OPTIONS}")
foreach(option IN LISTS SICONOS_CMAKE_OPTIONS)
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${option}")
endforeach()
if(WITH_COVERAGE)
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DWITH_TESTS_COVERAGE:BOOL=ON")
endif()
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")

if(SICONOS_COMPONENTS)
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

if(FROM_REPO)
  ctest_update()
endif()
ctest_configure(
  RETURN_VALUE CONFIGURE_RESULT
  CAPTURE_CMAKE_ERROR CONFIGURE_STATUS
  QUIET
  )

message("=============== End of ctest_configure =============== ")
message("------> Configure status/result : ${CONFIGURE_STATUS}/${CONFIGURE_RESULT}")

if(NOT CTEST_BUILD_CONFIGURATION)
  set(CTEST_BUILD_CONFIGURATION "Release")
endif()

message("\n\n=============== Start ctest_build =============== ")

ctest_build(
  PROJECT_NAME siconos-build
  CAPTURE_CMAKE_ERROR BUILD_STATUS
  RETURN_VALUE BUILD_RESULT
  QUIET
  )
message("=============== End of ctest_build =============== ")
message("------> Build status/result : ${BUILD_STATUS}/${BUILD_RESULT}")

# -- Tests --
message("\n\n=============== Start ctest_test (nbprocs = ${NP}) =============== ")
ctest_test(
  PARALLEL_LEVEL NP
  CAPTURE_CMAKE_ERROR TEST_STATUS
  SCHEDULE_RANDOM ON
  RETURN_VALUE TEST_RESULT
  QUIET
  )
message("=============== End of ctest_test =============== ")
message("------> Test status/result : ${TEST_STATUS}/${TEST_RESULT}")

ctest_test(PARALLEL_LEVEL ${NP} RETURN_VALUE TEST_RETURN_VAL)
message("=============== End of ctest_test =============== ")


if (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
  ctest_coverage()
endif (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  ctest_memcheck()
endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)

if(NOT SUBMIT)
  return()
endif()

# note: if the submission process experiences some slow-down, then we
# may get a return-code error, so we do it in a second phase.

# we need to get all the previously built files as ctest_start may
# begin with another tag
file(GLOB SUBMIT_FILES ${CMAKE_BINARY_DIR}/Testing/*/*)
message(STATUS "submit files : ${SUBMIT_FILES}")
ctest_submit(FILES ${SUBMIT_FILES}
  CAPTURE_CMAKE_ERROR  SUBMISSION_STATUS)

# ============= Summary =============
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "CTest process for siconos-build has ended.")
message(STATUS "Ctest model is: ${MODE}")
message(STATUS "Ctest executed on sources directory : ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY is: ${CTEST_BINARY_DIRECTORY}")
message(STATUS "CTEST_BUILD_CONFIGURATION is: ${CTEST_BUILD_CONFIGURATION}")
message(STATUS "Build name (cdash) : ${CTEST_BUILD_NAME}")
message(STATUS "Site (cdash) : ${CTEST_SITE}")
message(STATUS "=================================================================================================\n")

# tests failed?
if(NOT TEST_STATUS EQUAL 0 OR NOT TEST_RESULT EQUAL 0)
  message(FATAL_ERROR " *** test failure *** ")
endif()

# -- Submission failed? --
if(NOT SUBMISSION_STATUS EQUAL 0)
  message(WARNING " *** submission failure *** ")
endif()


if(NOT SUBMIT_RETURN_VAL EQUAL 0)
  message(WARNING " *** submission failure *** ")
endif()



