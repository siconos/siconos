# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.


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
    list(REMOVE_ITEM pyfiles occ_tools.py)
  endif()
  foreach(pyfile IN LISTS pyfiles)
    configure_file(${pyfile} ${SICONOS_PB11_BINARY_DIR}/siconos/io/${pyfile} @ONLY)
  endforeach()
  configure_file(format_v1.json ${SICONOS_PB11_BINARY_DIR}/siconos/io/format_v1.json COPYONLY)

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
