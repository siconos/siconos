# Cmake file used to for continuous integration to:
# - create docker container on Travis runner
# - configure, build and test Siconos on this runner.
#
# Called from python script (driver.py).

cmake_minimum_required(VERSION 3.7)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake;${CMAKE_SOURCE_DIR}/config)

include(Tools)

# Read and include siconos options from one or more files.
# 
if(CI_CONFIG)
  string(REPLACE "," ";" CI_CONFIG_LIST ${CI_CONFIG})
  foreach(_CI IN LISTS CI_CONFIG_LIST)
    include(config/${_CI}.cmake)
  endforeach()
endif()

message(STATUS "Siconos CMake options are : ${SICONOS_OPTIONS}")

# ---------- Docker or not docker ... ------------------------
if(WITH_DOCKER)
  find_program(DOCKER_COMMAND docker)
  if(DOCKER_COMMAND)
    set(DOCKER_MKSENV_COMMAND ${CMAKE_SOURCE_DIR}/machinery/mksenv.py)
    set(DOCKER_MKSENV_INPUT ${CMAKE_SOURCE_DIR}/config/siconos.yml)
    set(DOCKER_REPOSITORY ci_driver)
    # ctest driver file
    set(DOCKER_CTEST_DRIVER ${CMAKE_SOURCE_DIR}/ctest_driver_install_siconos.cmake)
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
