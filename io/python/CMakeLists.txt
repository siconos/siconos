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
  PARENT_NAME io
)

set_target_properties(pyio PROPERTIES OUTPUT_NAME "_io")

target_include_directories(pyio PRIVATE ${CMAKE_SOURCE_DIR}/io/src/mechanics)
target_include_directories(pyio PRIVATE ${CMAKE_SOURCE_DIR}/numerics/src/tools)

install(TARGETS
  pyio
  COMPONENT python
  LIBRARY DESTINATION "${SICONOS_PYTHON_INSTALL_DIR}/siconos/io"
)


if(HAVE_SICONOS_MECHANICS)
  # -- collect python files for the siconos/io module --
  collect_files(VAR pyfiles DIRS ${CMAKE_CURRENT_SOURCE_DIR} EXTS py)
  if(SICONOS_HAS_OpenCASCADE)
    message(STATUS "Siconos has OpenCASCADE")
  else()
    list(REMOVE_ITEM pyfiles occtools.py)
  endif()

  foreach(pyfile IN LISTS pyfiles)
    configure_file(${pyfile} ${SICONOS_PB11_BINARY_DIR}/siconos/io/${pyfile} @ONLY)
  endforeach()

  # -- collect python files used to create scripts --

  collect_files(VAR scriptfiles DIRS ${CMAKE_CURRENT_SOURCE_DIR}/scripts EXTS py)

  foreach(pyfile IN LISTS scriptfiles)
    get_filename_component(short_name ${pyfile} NAME_WE)
    configure_file(${pyfile} ${SICONOS_PB11_BINARY_DIR}/siconos/io/${short_name} @ONLY)
    install(PROGRAMS ${SICONOS_PB11_BINARY_DIR}/siconos/io/${short_name} DESTINATION bin)
  endforeach()

endif()



if(WITH_TESTING)
  #if(NOT WITH_BULLET)
  #  list(APPEND python_excluded_tests tests/test_mechanics_run.py)
  #endif()
  build_python_tests(EXCLUDE ${python_excluded_tests})
endif()
