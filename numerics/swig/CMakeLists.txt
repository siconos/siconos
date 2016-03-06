#================================================================
# cmake utilities to build python-bindings for numerics component
#================================================================

# wrap - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS
  ${SICONOS_SWIG_SRC_DIRS}/solverOptions.i
  ${SICONOS_SWIG_SRC_DIRS}/start.i
  ${SICONOS_SWIG_SRC_DIRS}/end.i
  ${SICONOS_SWIG_SRC_DIRS}/numpy.i)

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES "./${COMPONENT}")

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES)
# numerics is not written in C++
set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")
IF(NOT WITH_CXX)
  list(REMOVE_ITEM ${COMPONENT}_SWIG_DEFS "-dirvtable")
ENDIF(NOT WITH_CXX)

if(WITH_FCLIB)
  list(APPEND ${COMPONENT}_PYTHON_MODULES "./fclib")
else()
   list(APPEND ${COMPONENT}_python_excluded_tests
     ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_fclib.py)
endif()
include(swig_python_tools)

swig_module_setup(${COMPONENT}_PYTHON_MODULES)

build_python_tests()
