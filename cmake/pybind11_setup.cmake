# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2024 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#====================================================================================
# cmake utility to configure pybind11 (Python wrapper for C++) interface for Siconos. 
#====================================================================================
include(FindPythonModule)
message(STATUS "------------------------------------------------\n")
message(STATUS " Start pybind11 configuration to generate python packages for Siconos ...\n")
set(PYBIND11_FINDPYTHON ON) # Just to be sure. See https://pybind11.readthedocs.io/en/stable/cmake/index.html#modes
find_package(Python COMPONENTS Development Interpreter NumPy REQUIRED)

# Check if pybind11 is installed ...
find_python_module(pybind11 REQUIRED)

# And get its cmake root dir
execute_process(COMMAND 
  ${Python_EXECUTABLE} -m pybind11 --cmakedir # Get path to cmake-pybind11 config
  OUTPUT_VARIABLE pybind11_ROOT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

find_package(pybind11 CONFIG REQUIRED)

set(SICONOS_PB11_BINARY_DIR ${CMAKE_BINARY_DIR}/python
  CACHE INTERNAL "Working/binary for pybind11 and python stuff.")

if(WITH_TESTING)
  # Create test dir
  file(MAKE_DIRECTORY ${SICONOS_PB11_BINARY_DIR}/tests)
endif()

