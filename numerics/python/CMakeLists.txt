# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#
#


create_pybind11_module(
  MODULE_NAME pynumerics
  SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/numerics_pybind11.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lcp_pybind11.cpp ${CMAKE_CURRENT_SOURCE_DIR}/mlcp_pybind11.cpp ${CMAKE_CURRENT_SOURCE_DIR}/friction_contact_pybind11.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/relay_pybind11.cpp ${CMAKE_CURRENT_SOURCE_DIR}/generic_mechanical_pybind11.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/plasticity_pybind11.cpp
  DEPENDENCIES kernel numerics
  PARENT_NAME numerics
)

set_target_properties(pynumerics PROPERTIES OUTPUT_NAME "_numerics")

install(TARGETS
  pynumerics
  COMPONENT python
  LIBRARY DESTINATION "${SICONOS_PYTHON_INSTALL_DIR}/siconos/numerics"
)


collect_files(VAR pyfiles DIRS ${CMAKE_CURRENT_SOURCE_DIR} EXTS py RECURSIVE)

foreach(pyfile IN LISTS pyfiles)
  configure_file(${pyfile} ${SICONOS_PB11_BINARY_DIR}/siconos/numerics/${pyfile} @ONLY)
endforeach()


if(WITH_TESTING)
  list(APPEND python_excluded_tests tests/test_vi.py) # to be reviewed
  if(NOT SICONOS_HAS_MPI)
    list(APPEND python_excluded_tests tests/test_NM_MUMPS.py)
  endif()
  list(APPEND python_excluded_tests tests/test_fclib.py) # fclib pybind11 not yet implemented
  list(APPEND python_excluded_tests tests/test_lcp.py) # LCP pybind11 not yet implemented
  #list(APPEND python_excluded_tests tests/test_mc2d.py) # Now implemented!
  list(APPEND python_excluded_tests tests/test_mcp.py) # MCP pybind11 not yet implemented
  list(APPEND python_excluded_tests tests/test_mcp2.py) # MCP pybind11 not yet implemented
  list(APPEND python_excluded_tests tests/test_mlcp.py) # MLCP pybind11 not yet implemented
  list(APPEND python_excluded_tests tests/test_ncp.py) # NCP pybind11 not yet implemented
  list(APPEND python_excluded_tests tests/test_sparse.py) # TO BE REVIEWED/ pybind11 not fully implemented
  list(APPEND python_excluded_tests tests/test_tools.py) # Proj on cone pybind11 not yet implemented
  list(APPEND python_excluded_tests tests/test_linalg.py) # TO BE REVIEWED/ GFC3 pybind11 not fully implemented

  build_python_tests(EXCLUDE ${python_excluded_tests})
endif()
