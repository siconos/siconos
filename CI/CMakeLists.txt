# Configurations used for continuous integration
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake;${CMAKE_SOURCE_DIR}/config)

include(Tools)

# Check/set uninitialized cache variables
get_uninitialized_vars(CI_OPTIONS)

if(CI_CONFIG)
  string(REPLACE "," ";" CI_CONFIG_LIST ${CI_CONFIG})
  foreach(_CI ${CI_CONFIG_LIST})
    include(config/${_CI}.cmake)
  endforeach(_CI CI_CONFIG_LIST)
endif()

message(STATUS "CI options are : ${CI_OPTIONS}")
message(STATUS "Siconos CMake options are : ${SICONOS_OPTIONS}")

# ---------- Docker or not docker ... ------------------------
if(WITH_DOCKER)
  find_program(DOCKER_COMMAND docker)
  if(DOCKER_COMMAND)
    set(DOCKER_MKSENV_COMMAND ${CMAKE_SOURCE_DIR}/machinery/mksenv.py)
    set(DOCKER_MKSENV_INPUT ${CMAKE_SOURCE_DIR}/config/siconos.yml)
    set(DOCKER_REPOSITORY ci_driver)
    set(DOCKER_CTEST_DRIVER ${CMAKE_SOURCE_DIR}/machinery/CTestDriver.cmake)
    set(DOCKER_CTEST_MODE ${MODE})
    include(Docker)

    add_docker_targets(
      DOCKER_IMAGE ${DOCKER_IMAGE}
      DOCKER_TEMPLATE ${DOCKER_TEMPLATE}
      DOCKER_IMAGE_DIR ${DOCKER_IMAGE_DIR}
      DOCKER_REPOSITORY ${DOCKER_REPOSITORY}
      DOCKER_CTEST_DRIVER ${DOCKER_CTEST_DRIVER}
      DOCKER_SHARED_DIRECTORIES ${DOCKER_SHARED_DIRECTORIES}
      DOCKER_PROJECT_SOURCE_DIR ${DOCKER_PROJECT_SOURCE_DIR})

    project(ci_driver NONE)
  else()
    message(FATAL_ERROR "Cannot find docker.")
  endif()
endif()
