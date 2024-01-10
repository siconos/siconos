#-----------------------------------------------
# Ctest driver for siconos install.
# Target : continuous integration on gitlab-ci,
# aims at providing a proper install of siconos for a given configuration.
#
# Input variables :
# - USER_OPTIONS_FILE : user option file used by cmake to configure siconos. Default : siconos_conf.cmake.
#   Warning : always searched in siconos-tutorials/ci directory.
#   using common commands (lsb_release ...)
# ----------------------------------------------


# ============= setup  ================

# -- CI_PROJECT_DIR is a required environment variable --
# --> set by default for gitlab-ci, even inside the docker container

if(DEFINED ENV{GITLAB_CI})
  if($ENV{GITLAB_CI} STREQUAL true)
    set(CI_GITLAB ON)
  endif()
endif()
  
# if(NOT DEFINED ENV{CI_PROJECT_DIR} )
#   message(FATAL_ERROR "Please set env variable CI_PROJECT_DIR to siconos sources directory (git repo).")
# endif()


# - Source dir and path to siconos install
# if(NOT CTEST_SOURCE_DIRECTORY)
#   set(CTEST_SOURCE_DIRECTORY $ENV{CI_PROJECT_DIR})
# endif()

# # - Top level build directory -
# # If not specified : current dir.
# if(NOT CTEST_BINARY_DIRECTORY)
#   set(CTEST_BINARY_DIRECTORY .)
# endif()
# Build name (for cdash)


if(${CTEST_MODE} STREQUAL "configure" OR ${CTEST_MODE} STREQUAL "all")

  # -- Definition of all variables required for ctest --
  
  # Parallel build only for siconos_install. For examples it leads to: warning: jobserver unavailable: using -j1. Add `+' to parent make rule.
  #set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${CTEST_SOURCE_DIRECTORY}/cmake/valgrind.supp)
  
endif()

if(NOT CTEST_CMAKE_GENERATOR)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
endif()
  
include(${CTEST_SOURCE_DIRECTORY}/ci_gitlab/ctest_tools.cmake)

if(NOT CTEST_BUILD_NAME)
  set_cdash_build_name()
endif()

if(NOT CTEST_SITE)
  set_site_name()
endif()
  
# =============  Run ctest steps ================
# Either one by one (to split ci jobs) if CTEST_MODE=configure, build, test or
# all in a row if CTEST_MODE=all.
# Submit : only after test phase except if conf or build failed.

# - Configure -- 
if(${CTEST_MODE} STREQUAL "configure" OR ${CTEST_MODE} STREQUAL "all")

  # ---- CDash conf ----
  message("\n=============== Start conf for Siconos ctest/cdash ===============\n")
  # message("   - CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
  # message("   - CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")
  # message("   - CTEST_MODEL: ${model}")


  if(NOT CTEST_BUILD_CONFIGURATION)
    set(CTEST_BUILD_CONFIGURATION "Release")
  endif()

  # Write a note file for cdash server.
  # Content :
  # - info. regarding the runner, the system
  # - siconos config (user option file)
  write_notes() 

  # Current testing model. Priority: 
  # Nightly -> set by scheduler on gricad-gitlab
  # Continuous -> set in .gitlab-ci.yml
  # Experimental : default
  ctest_start(${model})

  # ----- Options to configure Siconos with cmake ---
  
  if(USER_OPTIONS_FILE)
    list(APPEND SICONOS_CMAKE_OPTIONS -DUSER_OPTIONS_FILE=${USER_OPTIONS_FILE})
  endif()

  list(APPEND SICONOS_CMAKE_OPTIONS -DWITH_GIT=ON) # required to generate siconos-commit.txt to tag cdash build in the examples.
  
  if(DEFINED ENV{OCE_INSTALL}) # set if oce has been installed using oce repo, in install_oce.sh
    message("Search oce in $ENV{OCE_INSTALL}.")
    list(APPEND SICONOS_CMAKE_OPTIONS -DOCE_DIR=$ENV{OCE_INSTALL})
  endif()

  # Set CTEST_CONFIGURE_COMMAND to cmake followed by siconos options 
  set(CTEST_CONFIGURE_COMMAND ${CMAKE_COMMAND})
  foreach(option IN LISTS SICONOS_CMAKE_OPTIONS)
    set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${option}")
  endforeach()
  set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")

  message("\n\n=============== Start ctest_configure =============== ")
  message("- Configure command line :\n ${CTEST_CONFIGURE_COMMAND}\n")

  ctest_configure(
    RETURN_VALUE _RESULT
    CAPTURE_CMAKE_ERROR _STATUS
    #QUIET
  )
  post_ctest(PHASE configure FORCE)
endif()
 
# - Build -
if(${CTEST_MODE} STREQUAL "build" OR ${CTEST_MODE} STREQUAL "all")
  
  if(${CTEST_MODE} STREQUAL "build")
    ctest_start(APPEND) # Restart from existing (configure step) cdash config
  endif()
  # --- Build ---

  message("\n\n=============== Start ctest_build =============== ")
  
  cmake_host_system_information(RESULT NP QUERY NUMBER_OF_LOGICAL_CORES)
  # if(NOT ALLOW_PARALLEL_BUILD)
  #   set(NP 1)
  # endif()
  ctest_build(
    FLAGS -j${NP}
    CAPTURE_CMAKE_ERROR _STATUS
    RETURN_VALUE _RESULT
  )
  post_ctest(PHASE build FORCE)
endif()

# - Test -
if(${CTEST_MODE} STREQUAL "test" OR ${CTEST_MODE} STREQUAL "all")
  # -- Tests --
  
  if(${CTEST_MODE} STREQUAL "test")
    ctest_start(APPEND)
  endif()
  message("\n\n=============== Start ctest_test =============== ")
  ctest_test(
    #PARALLEL_LEVEL NP
    CAPTURE_CMAKE_ERROR _STATUS
    #SCHEDULE_RANDOM ON
    RETURN_VALUE _RESULT
    # QUIET
  )
  post_ctest(PHASE test FORCE)

#   if(WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
#     #find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
#     # set(CTEST_COVERAGE_COMMAND gcov)

#     ctest_coverage(
#       CAPTURE_CMAKE_ERROR COVERAGE_STATUS
#       RETURN_VALUE COVERAGE_RESULT  
#       )
#   endif()

#   if(WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
#     ctest_memcheck()
#   endif()

# if(CDASH_SUBMIT)
#   ctest_submit(
#     RETURN_VALUE RETURN_STATUS
#     CAPTURE_CMAKE_ERROR SUBMISSION_STATUS
#   )
#   if(NOT SUBMISSION_STATUS EQUAL 0)
#     message(WARNING " *** submission failure *** ")
#   endif()
endif()

#   # -- memory check -- Skip this to 'enlight' submit process, since cdash inria is overbooked ...
#   # if(CTEST_BUILD_CONFIGURATION MATCHES "Profiling")
#   #   find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
#   #   set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-leak-kinds=definite,possible --track-origins=yes --error-limit=no --gen-suppressions=all") 
#   #   set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--quiet --leak-check=full --show-reachable=yes --error-limit=no --gen-suppressions=all") 
#   #   ctest_memcheck(PARALLEL_LEVEL NP QUIET)
#   # endif()

# endif()


# ============= Summary =============
message(STATUS "\n============================================ Summary ============================================")
message(STATUS "CTest process for Siconos has ended.")
message(STATUS "Ctest model is: ${model}")
message(STATUS "Ctest mode was: ${CTEST_MODE}")
message(STATUS "Ctest executed on sources directory : ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY is: ${CTEST_BINARY_DIRECTORY}")
message(STATUS "CTEST_BUILD_CONFIGURATION is: ${CTEST_BUILD_CONFIGURATION}")
if(CDASH_SUBMIT)
  message(STATUS "Cdash server name: ${CTEST_DROP_SITE}/${CTEST_DROP_LOCATION}.")
  message(STATUS "Cdash build name: ${CTEST_BUILD_NAME}")
  message(STATUS "Cdash site name: ${CTEST_SITE}")
endif()
message(STATUS "=================================================================================================\n")
