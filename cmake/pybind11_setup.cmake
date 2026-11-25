# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#
#====================================================================================
# cmake utility to configure pybind11 (Python wrapper for C++) interface for Siconos.
#====================================================================================
include(FindPythonModule)
message(" ---- Start pybind11 configuration to generate python packages for Siconos ... ----")
set(PYBIND11_FINDPYTHON ON) # Just to be sure. See https://pybind11.readthedocs.io/en/stable/cmake/index.html#modes
find_package(Python COMPONENTS Development Interpreter NumPy REQUIRED)

# Check if pybind11 is installed ...
# find_python_module(pybind11 REQUIRED)

# And get its cmake root dir
execute_process(COMMAND
  ${Python_EXECUTABLE} -m pybind11 --cmakedir # Get path to cmake-pybind11 config
  OUTPUT_VARIABLE pybind11_ROOT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

find_package(pybind11 3.0.2 CONFIG REQUIRED)
set(SICONOS_PB11_BINARY_DIR ${CMAKE_BINARY_DIR}/python
  CACHE INTERNAL "Working/binary for pybind11 and python stuff.")

if(WITH_TESTING)
  # Create test dir
  file(MAKE_DIRECTORY ${SICONOS_PB11_BINARY_DIR}/tests)
endif()
