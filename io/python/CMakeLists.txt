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
  MODULE_NAME pyio 
  SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/io_pybind11.cpp
  DEPENDENCIES kernel io
)

target_include_directories(pyio PRIVATE ${CMAKE_SOURCE_DIR}/io/src/mechanics)
target_include_directories(pyio PRIVATE ${CMAKE_SOURCE_DIR}/numerics/src/tools)

install(TARGETS
pyio
  COMPONENT python
  LIBRARY DESTINATION "${SICONOS_PYTHON_INSTALL_DIR}/siconos/io"
  )

# Installation des fichiers Python supplémentaires
install(FILES
    ${CMAKE_SOURCE_DIR}/io/swig/io/FrictionContactTrace.py
    ${CMAKE_SOURCE_DIR}/io/python/__init__.py
    DESTINATION "${SICONOS_PYTHON_INSTALL_DIR}/siconos/io")
