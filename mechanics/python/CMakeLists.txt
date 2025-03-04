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
#


create_pybind11_module(
  MODULE_NAME pymechanics 
  SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/mechanics_pybind11.cpp
  DEPENDENCIES kernel mechanics
)

target_include_directories(pymechanics PRIVATE ${CMAKE_SOURCE_DIR}/mechanics/src/collision/bullet)
target_include_directories(pymechanics PRIVATE ${CMAKE_SOURCE_DIR}/mechanics/src/collision/native/bodies)
target_include_directories(pymechanics PRIVATE ${CMAKE_SOURCE_DIR}/mechanics/src/collision/native)
target_include_directories(pymechanics PRIVATE ${CMAKE_SOURCE_DIR}/mechanics/src/collision)
target_include_directories(pymechanics PRIVATE ${CMAKE_SOURCE_DIR}/numerics/src/tools)

collect_files(VAR pyfiles DIRS ${CMAKE_CURRENT_SOURCE_DIR} EXTS py RECURSIVE)

foreach(pyfile IN LISTS pyfiles)
  configure_file(${pyfile} ${SICONOS_PB11_BINARY_DIR}/siconos/mechanics/${pyfile} @ONLY)
endforeach()


install(TARGETS
pymechanics
  COMPONENT python
  LIBRARY DESTINATION "${SICONOS_PYTHON_INSTALL_DIR}/siconos/mechanics"
  )
